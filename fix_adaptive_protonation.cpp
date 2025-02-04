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
enum {NEITHER = -1, SOLID = 0, SOLVENT = 1}

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::FixAdaptiveProtonation(LAMMPS* lmp, int narg, char** arg) : Fix(lmp, narg, arg), 
   pHstructureFile3(nullptr), pHstructureFile4(nullptr)
{
   if (narg < 10) utils::missing_cmd_args(FLERR, "fix adaptive_protonation", error);

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
      pHstructureFile3 = fopen(arg[4],"r"); // The structure in the solid state
      if (pHstructureFile1 == nullptr) error->one(FLERR,"Unable to open the file");
      pHstructureFile4 = fopen(arg[5],"r"); // The structure in the solvent phase ==> It should be modified based on the pKa and pH values by the fix constant_pH command
      if (pHstructureFile2 == nullptr) error->one(FLERR,"Unable to open the file");
   }

   typeOW = utils::numeric(FLERR,arg[6],false,lmp);

   threshold = utils::numeric(FLERR,arg[7],false,lmp);

   // Are these variables really necessary?
   pKa = utils::numeric(FLERR,arg[8],false,lmp);
   pH   = utils::numeric(FLERR,arg[9],false,lmp);
   req_numHs = utils::numeric(FLERR,arg[10],false,lmp);
   

}

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::~FixAdaptiveProtonation()
{
   if (pH1qs) memory->destroy(pH1qs);
   if (pH2qs) memory->destroy(pH2qs);
   if (typePerProtMol) memory->destroy(typePerProtMol);
   if (protonable) memory->destroy(protonable);
   if (mark) memory->destroy(mark);
   if (mark_prev) memory->destroy(mark_prev);
   if (mark_local) memory->destroy(mark_local);
   if (molecule_size) memory->destroy(molecule_size);
   if (molecule_size_local) memory->destroy(molecule_size_loca); 
   
   if (protonable_molids) delete [] protonable_molids;

   // It is a good practice to put the deallocated pointers to nullptr
   pH1qs = nullptr;
   pH2qs = nullptr;
   typePerProtMol = nullptr;
   protonable = nullptr;
   mark = nullptr;
   mark_prev = nullptr;
   mark_local = nullptr;
   molecule_size = nullptr;
   molecule_size_local = nullptr;

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
   nmolecules = atom->nmolecules;

   //It could be separate function
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

   if (atom->nmolecules > nmolecules)
   {
      memory->destroy(mark);
      memory->destroy(mark_prev);
      memory->destroy(mark_local);
      memory->destroy(molecule_size);
      memory->destroy(molecule_size_local);

      nmolecules = atom->nmolecules;
      
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

   mark_protonation_deprotonation();
   modify_protonable_hydrogens();
   set_protonable_molids();
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
       pHnStructures2 = s9es2,1,MPI_INT,0,world);
   }

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
   Setting the molids for protonable molecules
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::set_protonable_molids()
{
   int nlocal = atom->nlocal;
   int * molecule = atom->molecule;
   int n_protonable_local = 0;
   n_protonable = 0;
   int * protonable_molids_local = new int[nlocal];
   int * protonable_molids;

   int nprocs = comm->nprocs;
   int * recv_counts = new int[nprocs];
   int * displs = new int[nprocs];

   
   for (int i = 0; i < nlocal; i++) {
      if (mark[i] == 1) {
         protonable_molids_local[n_protonable_local++] = molecule[i];
      }
   }


   /* Extracting the unique molecule_ids so that I would know the 
      n_molecule_ids and the required size of the molids array*/
   
    // Collect all molecule IDs from all processes
    MPI_Allreduce(&n_protonable_local,&n_protonable, 1, MPI_INT, MPI_SUM, world);

    protonable_molids = new int[n_protonable];
    

    MPI_Allgather(&n_protonable, 1, MPI_INT, recv_counts, 1, MPI_INT, world);

    displs[0] = 0;
    for (int i = 1; i < nprocs; i++) {
        displs[i] = displs[i - 1] + recv_counts[i - 1];
    }

    MPI_Allgatherv(protonable_molids_local, n_protonable_local, MPI_INT, protonable_molids, recv_counts, displs, MPI_INT, world);

    std::set<int> unique_molids(protonable_molids, protonable_molids + n_protonable);

    delete [] protonable_molids;

    n_protonable = unique_molids.size();

    protonable_molids = new int[n_protonable];

    int i = 0;
    for (const int &mol_id : unique_molids) {
        protonable_molids[i++] = mol_id;
    }

    delete [] protonable_molids_local;
    delete [] recv_counts;
    delete [] displs;
   
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


   // You have take care of the situation in which one of the atoms of a molecule has mark[i] == 1
   for (int i; i < nlocal; i++) {
      switch(mark[molecule[i]])
      {
         case NEITHER: // Not protonable --> nothing to do here
            break;
         case SOLVENT: // In the water
            if (mark_prev[molecule[i]] == SOLID) // It came from the solid --> protonate it
               q[i] = pH2qs[type[i]][0];
            else
            {} // It used to be in the water --> do nothing
            break;
         case SOLID: // In the solid
            if (make_prev[molecule[i]] == SOLVENT) // It came from the water --> deprotonate it
               q[i] = pH1qs[type[i]][0];
            else
            {} // It used to be in the solid --> do nothing 
      }
   }
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
   return static_cast<double>(natoms_change_total);
}

/* --------------------------------------------------------------------------
   Output the changes in the topology --> nbonds and nangles
   -------------------------------------------------------------------------- */

double FixAdaptiveProtonation::compute_vector(int n)
{
   switch (n)
   {
      case 1:
         return static_cast<double>(nbonds_change_total);
      case 2:
         return static_cast<double>(nangles_change_total);
   }
   error->all(FLERR,"Undefined error");
}

/* --------------------------------------------------------------------------
   This part needs to be updated in the final version ....
   -------------------------------------------------------------------------- */
   
double FixAdaptiveProtonation::memory_usage()
{
   return 0.0;
}


