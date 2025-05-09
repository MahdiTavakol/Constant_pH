/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ---------- v0.10.00----------------- */
// Please remove unnecessary includes 
#include "fix_adaptive_protonation.h"

#include "atom.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include "region.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#include "respa.h"
#include "update.h"
#include "variable.h"


#include "force.h"
#include "pair.h"
#include "improper.h"
#include "dihedral.h"
#include "kspace.h"
#include "angle.h"
#include "bond.h"
#include "atom.h"
#include "group.h"

#include "thermo.h"
#include <cmath>
#include <cstring>
#include <stdio.h>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { NONE, CONSTANT, EQUAL, ATOM };
enum {NEITHER = -1, SOLID = 0, SOLVENT = 1};

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::FixAdaptiveProtonation(LAMMPS* lmp, int narg, char** arg) : Fix(lmp, narg, arg), 
   pHStructureFile1(nullptr), pHStructureFile2(nullptr)
{
   if (narg < 7) utils::missing_cmd_args(FLERR, "fix adaptive_protonation", error);

   dynamic_group_allow = 0;
   scalar_flag = 1; 
   vector_flag = 1;
   size_vector = 2;
   global_freq = 1; // What is this?
   extscalar = 1; // What is this?
   extvector = 1; // What is this?
   energy_global_flag = 1;
   virial_global_flag = 1;
      int typeP, typeO, typeH, typeOH, typePOH;
   virial_peratom_flag = 1;
   respa_level_support = 1;
   
   nevery = utils::numeric(FLERR, arg[3], false, lmp);

   if (nevery < 0) error->all(FLERR,"Illegal fix adaptive_protonation every value {}", nevery);

   if (comm->me == 0) {
      pHStructureFile1 = fopen(arg[4],"r"); // The structure in the solid state
      if (pHStructureFile1 == nullptr) error->one(FLERR,"Unable to open the file");
      pHStructureFile2 = fopen(arg[5],"r"); // The structure in the solvent phase ==> It should be modified based on the pKa and pH values by the fix constant_pH command
      if (pHStructureFile2 == nullptr) error->one(FLERR,"Unable to open the file");
   }

   typeOW = utils::numeric(FLERR,arg[6],false,lmp);
   threshold = utils::numeric(FLERR,arg[7],false,lmp);

   molecule_id_flag = false;
   int iarg = 8;
   while (iarg < narg) {
      if (strcmp(arg[iarg], "reset_molecule_ids") == 0) {
         molecule_id_flag = true;
         iarg++;
      }
      else
         error->all(FLERR,"Unknown keyword");
   }

}

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::~FixAdaptiveProtonation()
{
   if (pH1qs) memory->destroy(pH1qs);
   if (pH2qs) memory->destroy(pH2qs);
   if (typePerProtMol) memory->destroy(typePerProtMol);
   if (protonable) memory->destroy(protonable);
   
   if (protonable_molids) delete [] protonable_molids;

   deallocate_storage(); 

   // It is a good practice to put the deallocated pointers to nullptr
   pH1qs = nullptr;
   pH2qs = nullptr;
   typePerProtMol = nullptr;
   protonable = nullptr;

   protonable_molids = nullptr;
}

/* --------------------------------------------------------------------------------------- */

int FixAdaptiveProtonation::setmask()
{
   int mask = 0;
   mask |= PRE_EXCHANGE;
   return mask;
}

/* --------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::init()
{
   // I am not sure if this should be here to in the setup() function
   read_pH_structure_files();

   // Checking if the atom style contains the molecules information
   if (atom->molecular != 1) error->all(FLERR,"Illegal atom style in the fix adpative protonation");

   nmax = atom->nmax;
   nmolecules = atom->nmolecule;

   allocate_storage();

   if (molecule_id_flag)
      set_molecule_id();
}

/* ---------------------------------------------------------------------------------------
    It is need to access the neighbor list
   --------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::init_list(int /*id*/, NeighList* ptr) 
{
   list = ptr;
}

/* --------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::pre_exchange()
{
   if(update->ntimestep != next_reneighbor) return;

   if (atom->nmolecule > nmolecules)
   {
      nmolecules = atom->nmolecule;
      deallocate_storage();
      allocate_storage();
   }

   // Counting the number of water molecules surrounding the protonable molecules
   mark_protonation_deprotonation();

   // This is required since the fix_constant_pH.cpp does not deal with those molecules in the solid
   modify_protonation_state();

   // Resetting the mark_prev parameter to help us keep the track of which molecule moves from solid to solvent and vice versa
   set_mark_prev();

}

/* ----------------------------------------------------------------------------------------
    Reading the structure files
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::read_pH_structure_files()
{
    /* File format
    * Comment
    * pHnStructures
    * pHnTypes
    * type1,  number of type1 atoms in the protonable molecule, qState1, qState2, qState3
    * ... 
    * ...  
    */

   /*Allocating the required memory*/
   int ntypes = atom->ntypes;
   memory->create(protonable,ntypes+1,"constant_pH:protonable"); //ntypes+1 so the atom types start from 1.
   memory->create(typePerProtMol,ntypes+1,"constant_pH:typePerProtMol");



   char line[128];
   if (comm->me == 0)
   {
       if (!pHStructureFile1 || !pHStructureFile2 )
           error->all(FLERR,"Error in reading the pH structure file in fix constant_pH");

       // comment 
       fgets(line,sizeof(line),pHStructureFile1);
       fgets(line,sizeof(line),pHStructureFile2);

       // pHnStructures
       fgets(line,sizeof(line),pHStructureFile1);
       line[strcspn(line,"\n")] = '\0';
       char *token = strtok(line,",");   
       pHnStructures1 = std::stoi(token);

       // pHnStructures
       fgets(line,sizeof(line),pHStructureFile2);
       line[strcspn(line,"\n")] = '\0';
       token = strtok(line,",");   
       pHnStructures2 = std::stoi(token);
   }

   MPI_Bcast(&pHnStructures1,1,MPI_INT,0,world);
   MPI_Bcast(&pHnStructures2,1,MPI_INT,0,world);

   memory->create(pH1qs,ntypes+1,pHnStructures1, "constant_pH:pH1qs");
   memory->create(pH2qs,ntypes+1,pHnStructures2, "constant_pH:pH2qs");

   if (comm->me == 0) {
      // pHnTypes
      fgets(line,sizeof(line),pHStructureFile1);
      line[strcspn(line,"\n")] = '\0';
      char *token = strtok(line,",");    
      pHnTypes1 = std::stoi(token);

      fgets(line,sizeof(line),pHStructureFile2);
      line[strcspn(line,"\n")] = '\0';
      token = strtok(line,",");    
      pHnTypes2 = std::stoi(token);
	
      for (int i = 1; i < ntypes+1; i++)
      {
	        protonable[i] = 0;
	        typePerProtMol[i] = 0;
	        for (int j = 0; j < pHnStructures1; j++)
	           pH1qs[i][j] = 0.0;
	        for (int j = 0; j < pHnStructures2; j++)
	           pH2qs[i][j] = 0.0;
      }  
	   
      for (int i = 0; i < pHnTypes1; i++)
      {
	        if (fgets(line,sizeof(line),pHStructureFile1) == nullptr)
	           error->all(FLERR,"Error in reading the pH structure file in fix constant_pH");
	        line[strcspn(line,"\n")] = '\0';
	        token = strtok(line,",");
	        int type = std::stoi(token);
	        protonable[type] = 1;
	        token = strtok(NULL,",");
	        typePerProtMol[type] = std::stoi(token);
	        for (int j = 0; j < pHnStructures1; j++) {
	           token = strtok(NULL,",");
	           pH1qs[type][j] = std::stod(token);  
	        }

	        if (fgets(line,sizeof(line),pHStructureFile2) == nullptr)
	           error->all(FLERR,"Error in reading the pH structure file in fix constant_pH");
	        line[strcspn(line,"\n")] = '\0';
	        token = strtok(line,",");
	        type = std::stoi(token);
	        protonable[type] = 1;
	        token = strtok(NULL,",");
	        typePerProtMol[type] = std::stoi(token);
	        for (int j = 0; j < pHnStructures2; j++) {
	           token = strtok(NULL,",");
	           pH2qs[type][j] = std::stod(token);  
	         }
       }
       fclose(pHStructureFile1);
       fclose(pHStructureFile2);
   }
   
   // To avoid the code trying to reclose it.
   pHStructureFile1 = nullptr;
   pHStructureFile2 = nullptr;
   
   MPI_Bcast(protonable,ntypes+1,MPI_INT,0,world);
   MPI_Bcast(typePerProtMol,ntypes+1,MPI_INT,0,world);
   MPI_Bcast(pH1qs[0],(ntypes+1)*(pHnStructures1),MPI_DOUBLE,0,world);
   MPI_Bcast(pH2qs[0],(ntypes+1)*(pHnStructures2),MPI_DOUBLE,0,world); 
}

/* ----------------------------------------------------------------------------------------
   Deallocating the storage

   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::deallocate_storage()
{
   memory->destroy(mark);
   memory->destroy(mark_prev);
   memory->destroy(mark_local);
   memory->destroy(molecule_size);
   memory->destroy(molecule_size_local);
   mark = nullptr;
   mark_prev = nullptr;
   mark_local = nullptr;
   molecule_size = nullptr;
   molecule_size_local = nullptr;
}

/* ----------------------------------------------------------------------------------------
   Allocating the storage

   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::allocate_storage()
{
   memory->create(mark,nmolecules+1,"AdaptiveProtontation:mark");
   memory->create(mark_prev,nmolecules+1,"AdaptiveProtonation:mark_prev");
   memory->create(mark_local,nmolecules+1,"AdaptiveProtonation:mark_local");
   memory->create(molecule_size,nmolecules+1,"AdaptiveProtonation:molecule_size");
   memory->create(molecule_size_local,nmolecules+1,"AdaptiveProtonation:molecule_size_local");
   std::fill(mark,mark+nmolecules+1,0);
   std::fill(mark_prev,mark_prev+nmolecules+1,-1); // I put it on purpose so in the first step every molecule changes
   std::fill(mark_local,mark_local+nmolecules+1,0);
   std::fill(molecule_size,molecule_size+nmolecules+1,0);
   std::fill(molecule_size_local,molecule_size_local+nmolecules+1,0);
}

/* ----------------------------------------------------------------------------------------
   getting the number of water molecules near phosphates and 
   flag them for protonation/deprotonation
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::mark_protonation_deprotonation()
{
   int *ilist, *jlist, *numneigh, **firstneigh;
   int inum, jnum;
   int wnum; // number of surrounding water molecules

   inum = list->inum; // I do not ghost atoms for inum. however, I need them in jnum
   ilist = list->ilist;
   numneigh = list->numneigh;
   firstneigh = list->firstneigh;

   int * type = atom->type;
   int * molecule = atom->molecule;

   for (int ii = 0; ii < inum; ii++) {
      wnum = 0.0;
      int i = ilist[ii];
      molecule_size_local[molecule[i]]++;

      // Check if this atom is protonable --> if not do not bother with it.
      if (protonable[type[i]] == 0) {
         mark_local[molecule[i]] = NEITHER;
         continue;
      }
      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++) {
         int j = jlist[jj];
         j &= NEIGHMASK;

         if (type[j] == typeOW)
            wnum++;  // Just considering the Oxygens. It is possible that both O and H from the same water molecule are close to this atom.
      }
      if (wnum >= threshold) {
         mark_local[molecule[i]] += SOLVENT;
      } else {
         mark_local[molecule[i]] += SOLID;
      }
   }


   // Reducing the values from various cpus
   MPI_Allreduce(&mark_local,&mark,nmolecules+1,MPI_DOUBLE,MPI_SUM,world);
   MPI_Allreduce(&molecule_size_local,&molecule_size,nmolecules+1,MPI_DOUBLE,MPI_SUM,world);

   for (int i = 1; i < nmolecules+1; i++)
   {
      double test_condition = mark[i]/molecule_size[i];
      if (test_condition >= 0 && test_condition <= 0.5 ) mark[i] = SOLID;
      else if (test_condition > 0.5 && test_condition <= 1) mark[i] = SOLVENT;
      else if (test_condition > 1 || test_condition < -1) 
         error->one(FLERR,"Error in fix adaptive_protonation: You should never have reached here!");
      else mark[i] = NEITHER;
   }  
     
}

/* ----------------------------------------------------------------------------------------
   Setting separate molecule ids for different phosphate ions 
   It might need to be a separate command in LAMMPS
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::set_molecule_id()
{
   int natom = atom->nlocal + atom->nghost;
   int nlocal = atom->nlocal;
   int nmax = atom->nmax;
   int * molecule = atom->molecule;
   int * num_bond = atom->num_bond;
   int ** bond_atom = atom->bond_atom;
   int * tag = atom->tag; // atom-id



   for (int i = 0; i < nlocal; i++) {
      for (int k = 0; k < num_bond[i]; k++) {
         int jtag = bond_atom[i][k]; // the tag (atom-id) of kth bonds of atom i
         int j = atom->map(jtag);
         if (j == -1) {
            error->warning(FLERR,"Bond atom missing in fix AdaptiveProtonation");
            continue;
         }
         molecule[i] = MIN(molecule[i],molecule[j]); // I am not sure about header for the MIN
         molecule[j] = molecule[i];
      }
   }

   // You need to think about neighbor exchange;
   /*
   no need for an exchange 
   as LAMMPS itself takes care of exchange.

   if atom i is in the proc n connected to atom j in proc n both of the atoms 
   will be among the ghost atoms of the other atom. And as the minimum of the molecule_id
   is the same in both the procs both the atoms would end up having the same 
   molecule id and so there is no need for an atom exchange.
   Of course if the molecules are long spanning multiple procs there is a need
   for atom exchange here.
   */
}

/* ----------------------------------------------------------------------------------------
   Getting the molids for protonable molecules
   The molids must be allocated otherwise an error occurs
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::get_protonable_molids(int *_molids) const 
{
   if (_molids == nullptr) error->all(FLERR,"The _molids array in the fix adaptive protonation get protonable_molids must be allocated");

   for (int i = 0; i < n_protonable; i++) {
      _molids[i] = protonable_molids[i];
   }
}

/* ----------------------------------------------------------------------------------------
   Changing from the protonated to deprotonated states --> Moving from the solvent to the solid phase

   ---------------------------------------------------------------------------------------- */
void FixAdaptiveProtonation::modify_protonation_state()
{
   int nlocal = atom->nlocal;
   double* q = atom->q;
   int * type = atom->type;
   int * molecule = atom->molecule;
   int nchanges_local[3] = {0,0,0};


   // You have take care of the situation in which one of the atoms of a molecule has mark[i] == 1
   for (int i; i < nlocal; i++) {
      switch(mark[molecule[i]])
      {
         case NEITHER: // Not protonable --> nothing to do here
            break;
         case SOLVENT: // In the water
            // It came from the solid --> protonate it
            if (mark_prev[molecule[i]] == SOLID) {
               q[i] = pH2qs[type[i]][0]; // I just chose the first structure the fix constant pH itself selects the appropriate one
               nchanges_local[0]++;
               nchanges_local[1]++;
            }
            else
            {} // It used to be in the water --> do nothing
            break;
         case SOLID: // In the solid
            // It came from the water --> deprotonate it
            if (mark_prev[molecule[i]] == SOLVENT) {
               q[i] = pH1qs[type[i]][0]; // I would expect that there is just one structure in the solid state
               nchanges_local[0]++;
               nchanges_local[2]++;
            }
            else
            {} // It used to be in the solid --> do nothing 
      }
   }

   MPI_Allreduce(nchanges_local,nchanges,3,MPI_INT,MPI_SUM,world);
}


/* --------------------------------------------------------------------------
   Set the mark_prev
   -------------------------------------------------------------------------- */

void FixAdaptiveProtonation::set_mark_prev()
{
   for (int i = 0; i < nmolecules + 1; i++)
      mark_prev[i] = mark[i];
}

/* --------------------------------------------------------------------------
   Output the changes in the number of hydrogen atoms
   -------------------------------------------------------------------------- */

double FixAdaptiveProtonation::compute_scalar()
{
}

/* --------------------------------------------------------------------------
   Output the changes in the topology --> nbonds and nangles
   -------------------------------------------------------------------------- */

double FixAdaptiveProtonation::compute_vector(int n)
{
   if (n < 3)
      return nchanges[n];
   else
      error->one(FLERR,"Out of bound access");
}

/* --------------------------------------------------------------------------
   This part needs to be updated in the final version ....
   -------------------------------------------------------------------------- */
   
double FixAdaptiveProtonation::memory_usage()
{
   return 0.0;
}


