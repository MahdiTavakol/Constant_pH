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
/* ---------- v0.05.07----------------- */
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

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::FixAdaptiveProtonation(LAMMPS* lmp, int narg, char** arg) :
   Fix(lmp, narg, arg), typePstr(nullptr), typeOstr(nullptr), typeHstr(nullptr),
   typeOHstr(nullptr), typePOHstr(nullptr)
{
   if (narg < 11) utils::missing_cmd_args(FLERR, "fix AdaptiveProtonation", error);

   dynamic_group_allow = 0;
   scalar_flag = 1; 
   vector_flag = 1;
   size_vector = 2;
   global_freq = 1; // What is this?
   extscalar = 1; // What is this?
   extvector = 1; // What is this?
   energy_global_flag = 1;
   virial_global_flag = 1;
   virial_peratom_flag = 1;
   respa_level_support = 1;
   
   nevery = utils::numeric(FLERR, arg[3], false, lmp);

   if (utils::strmatch(arg[4], "^v_")) {
      typePstr = utils::strdup(arg[4]+2);
   }
   else
   {
      typeP = utils::numeric(FLERR,arg[4],false,lmp);
      typePstyle = CONSTANT;
   }
   if (utils::strmatch(arg[5], "^v_")) {
      typeOstr = utils::strdup(arg[5]+2);
   }
   else
   {
      typeO = utils::numeric(FLERR,arg[5],false,lmp);
      typeOstyle = CONSTANT;
   }
   if (utils::strmatch(arg[6], "^v_")) {
      typeHstr = utils::strdup(arg[6]+2);
   }
   else
   {
      typeH = utils::numeric(FLERR,arg[6],false,lmp);
      typeHstyle = CONSTANT;
   }
   if (utils::strmatch(arg[7], "^v_")) {
      typeOHstr = utils::strdup(arg[7]+2);
   }
   else
   {
      typeOH = utils::numeric(FLERR,arg[7],false,lmp);
      typeOHstyle = CONSTANT;
   }
   if (utils::strmatch(arg[8], "^v_")) {
      typePOHstr = utils::strdup(arg[8]+2);
   }
   else {
      typePOH = utils::numeric(FLERR,arg[8],false,lmp);
      typePOHstyle = CONSTANT;  
   }
   if (utils::strmatch(arg[8], "^v_")) {
      error->all(FLERR,"The variable input for threshold in AdaptiveProtonation is not supported yet!");
   }
   else
   {
      threshold = utils::numeric(FLERR,arg[8],false,lmp);
   }
   pKa = utils::numeric(FLERR,arg[9],false,lmp);
   pH   = utils::numeric(FLERR,arg[10],false,lmp);
   req_numHs = utils::numeric(FLERR,arg[11],false,lmp);
   
   nmax = 1;
   memory->create(mark,nmax,"AdaptiveProtontation:mark");
}

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::~FixAdaptiveProtonation()
{
   if (typePstr) delete[] typePstr;
   if (typeOstr) delete[] typeOstr;
   if (typeHstr) delete[] typeHstr;
   if (typeOHstr) delete[] typeOHstr;
   if (typePOHstr) delete[] typePOHstr;
   if (protonable_molids) delete[] protonable_molids;

   memory->destroy(mark);
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
   // check variables
   if (typePstr) {
      typePvar = input->variable->find(typePstr);
      if (typePvar < 0) error->all(FLERR, "Variable {} for fix AdaptiveProtonation does not exist", typePstr);
      if (input->variable->equalstyle(typePvar))
          typePstyle = EQUAL;
      else if (input->variable->atomstyle(typePvar))
          error->all(FLERR, "Atomic style variable is not supported in fix adaptiveProtonaton");
      else
          error->all(FLERR, "Variable {} for fix adaptiveProtonatonis invalid style", typePvar);
   }
   if (typeOstr) {
      typeOvar = input->variable->find(typeOstr);
      if (typeOvar < 0) error->all(FLERR, "Variable {} for fix AdaptiveProtonation does not exist", typeOstr);
      if (input->variable->equalstyle(typeOvar))
          typeOstyle = EQUAL;
      else if (input->variable->atomstyle(typeOvar))
          error->all(FLERR, "Atomic style variable is not supported in fix adaptiveProtonaton");
      else
          error->all(FLERR, "Variable {} for fix adaptiveProtonatonis invalid style", typeOvar);
   }
   if (typeHstr) {
      typeHvar = input->variable->find(typeHstr);
      if (typeHvar < 0) error->all(FLERR, "Variable {} for fix AdaptiveProtonation does not exist", typeHstr);
      if (input->variable->equalstyle(typeHvar))
          typeHstyle = EQUAL;
      else if (input->variable->atomstyle(typeHvar))
          error->all(FLERR, "Atomic style variable is not supported in fix adaptiveProtonaton");
      else
          error->all(FLERR, "Variable {} for fix adaptiveProtonatonis invalid style", typeHvar);
   }
   if (typeOHstr) {
      typeOHvar = input->variable->find(typeOHstr);
      if (typeOHvar < 0) error->all(FLERR, "Variable {} for fix AdaptiveProtonation does not exist", typeOHstr);
      if (input->variable->equalstyle(typeOHvar))
          typeOHstyle = EQUAL;
      else if (input->variable->atomstyle(typeOHvar))
          error->all(FLERR, "Atomic style variable is not supported in fix adaptiveProtonaton");
      else
          error->all(FLERR, "Variable {} for fix adaptiveProtonatonis invalid style", typeOHvar);
   }
   if (typePOHstr) {
      typePOHvar = input->variable->find(typePOHstr);
      if (typePOHvar < 0) error->all(FLERR,"Variable {} for fix AdaptiveProtonation does not exit", typePOHstr);
      if (input->variable->equalstyle(typePOHvar))
         typePOHstyle = EQUAL;
      else if (input->variable->atomstyle(typePOHvar))
         error->all(FLERR,"Atomic style variable is not supported in fix AdaptiveProtonation");
      else
         error->all(FLERR,"Variable {} for fix adaptiveProtonation is invalid style",typePOHvar);
   }


   // reseting the molecule ids so each atom and all of its connected atoms have the same molecule id
   int * molecule = atom->molecule;
   int * tag = atom->tag;
   int nlocal = atom->nlocal;
   for (int i = 0; i < nlocal; i++) {
      molecule[i] = tag[i];
   }
   set_molecule_id();

   // The initial value for the change in the natoms, nbonds and nangles
   natoms_change = 0;
   nbonds_change = 0;
   nangles_change = 0;
   natoms_change_total = 0;
   nbonds_change_total = 0;
   nangles_change_total = 0;
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

   if (atom->nmax > nmax)
   {
      memory->destroy(mark);
      nmax = atom->nmax;
      memory->create(mark,nmax,"AdaptiveProtontation:mark");
   }

   mark_protonation_deprotonation();
   modify_protonable_hydrogens();
   set_protonable_molids();
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

   inum = list->inum + list->gnum;
   ilist = list->ilist;
   numneigh = list->numneigh;
   firstneigh = list->firstneigh;

   int * type = atom->type;

   for (int ii = 0; ii < inum; ii++) {
      wnum = 0.0;
      int i = ilist[ii];
      if (type[i] != typeP) {
         mark[i] = 0.0;
         continue;
      }
      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++) {
         int j = jlist[jj];
         j &= NEIGHMASK;

         if (type[j] == typeO)
            wnum++;  // Just considering the Oxygens. It is possible that both O and H from the same water molecule are close to this atom.
      }
      if (wnum >= threshold) {
         mark[i] = 1;
      }
      else {
         mark[i] = -1;
      }  
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
   changing the protonation state of phosphates 
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::modify_protonable_hydrogens()
{
   int nlocal = atom->nlocal;
   tagint *tag = atom->tag;
   int * type = atom->type;
   double **x = atom->x;
   AtomVec *avec = atom->avec;


   if (atom->avec->bonds_allow) {
      
      int pAtom;     // local atom id for the P atom
      int oAtoms[3]; // local atom id for the O atoms  <----> It is only used for identifying those oxygen atoms needed to be added
      int hAtoms[3]; // local atom id for the H atoms  <----> It is only used for identifying those hydrogen atoms needed to be removed

      for (int i = 0; i < nlocal; i++) {
         if (mark[i] == 0) continue;


         pAtom = i;
         int numHs = 0; 
         for (int m = 0; m < atom->num_bond[i]; m++) {
            /* 
               bond_atoms[i][m] contains the global id (tag) of the connected atom
               so there is a need for mapping the global id to the local one
            */
            int atom2 = atom->map(atom->bond_atom[i][m]);
            if (m < 3) oAtoms[m] = atom2;
            else error->warning(FLERR,"Number of atoms connected to P is {}",atom->num_bond[i]);
            if (atom2 == -1) error->warning(FLERR,"Bond atom missing in fix AdaptiveProtonation");
            for (int n = 0; n < atom->num_bond[atom2]; n++) {
               int atom3 = atom->map(atom->bond_atom[atom2][n]);
               if (type[atom3] == typeH) 
               {
                   if (numHs < 3)
                      hAtoms[numHs++] = atom3;
                   else
                      error->warning(FLERR,"Number of protons are {}",++numHs);
               }
            }

            /* 
               Considering different values of mark[i] and setting the number of hydrogen atoms
               to an appropriate value
            */
            
            if (mark[i] == -1 && numHs == 0) continue; // PO4 inside the HAp ---> Nothing to do here!
            if (mark[i] == -1 && numHs > 0) remove_hydrogens(i,numHs,oAtoms,hAtoms); //HPO4, H2PO4 or H3PO4 inside ---> The hydrogens must be removed
            if (mark[i] == 1 && numHs == req_numHs) continue; // HxPO4 inside the surface with the exact number of required Hydrogen atoms ---> Nothing to do here!
            if (mark[i] == 1 && numHs > req_numHs) remove_hydrogens(i,numHs-req_numHs,oAtoms,hAtoms); // HxPO4 in the surface or in solution with higher number of Hydrogens
            if (mark[i] == 1 && numHs < req_numHs) {
               /* HxPO4 in the surface or in solution with lower number of Hydrogens
                  I would prefer to remove all the hydrogens of this phosphate and then add
                  the required number to eliminate the possibility of having an Oxygen atom with
                  two hydrogen atoms
               */
               remove_hydrogens(i,numHs,oAtoms, hAtoms);
               add_hydrogens(i,req_numHs,oAtoms, hAtoms);
            }
         }
      }
   }
}

/* ------------------------------------------------------------------------------------------
   Adding the requested number of hydrogen atoms to the molecule with P local id of i
   Warning <----> This function does not take into account the number of hydrogen atoms
   aldready present in the phosphate molecule so the best practice is to remove all
   the hydrogen atoms and then add the required number of hydrogens
   ------------------------------------------------------------------------------------------ */

void FixAdaptiveProtonation::add_hydrogens(const int& i, const int& req_numHs, const int oAtoms[3], const int hAtoms[3]) 
{
   double **x = atom->x;

   // adding hydrogens to this phosphate
   double coors[3][3]; // Maximum three phosphate ions
               
   coors[0][0] = x[oAtoms[0]][0] + 0.9;
   coors[0][1] = x[oAtoms[0]][1];
   coors[0][2] = x[oAtoms[0]][2];
   coors[1][0] = x[oAtoms[1]][0] - 0.8;
   coors[1][1] = x[oAtoms[1]][1];
   coors[1][2] = x[oAtoms[1]][2];
   coors[2][0] = x[oAtoms[2]][0];
   coors[2][1] = x[oAtoms[2]][1] + 0.64;
   coors[2][2] = x[oAtoms[2]][2] - 0.64;
               

   // atom->avec->create_atoms itself updates the atom->nlocal value

               
   // Let's connect these two new Hs to the P
   if (atom->num_bond[i] == atom->bond_per_atom)
      error->one(FLERR,"Num bonds exceeded bonds per atom in fix AdaptiveProtonation");

   /*
     I guess that the size of bond_type and bond_atom is atom->bond_per_atom so there is enough space there
   */

   for (int k = 0; k < req_numHs; k++) {
      // Adding the kth hydrogen atom
      
      error->warning(FLERR,"atom->natoms before adding atoms = {}",atom->natoms);
      atom->avec->create_atom(typeH,coors[k]);
      error->warning(FLERR,"atom->natoms after adding atoms = {}",atom->natoms);
      // If these two are the same I need to modify the atom->nlocal ---> However, it is not so likely.
      natoms_change++;
      

      // global and local indexes of new hydrogens are needed to modify the bond and angle information in the topology
      int tagH = (atom->tag[atom->nlocal-1] > 0)? atom->tag[atom->nlocal-1]: atom->natoms-1;
      int jH = atom->nlocal-1 ;

      
      // The kth OH bond and angle (Fortunately each hydrogen contributes to just one bond and one angle!)
      int p_local_index = i;
      int o_local_index = oAtoms[k];
      int h_local_index = jH;
      int p_global_index = atom->tag[p_local_index];
      int o_global_index = atom->tag[o_local_index];
      int h_global_index = tagH;
      int & p_num_bond = atom->num_bond[p_local_index]; // it is here just for completeness, we are not adding any new bonds to the P atom here.
      int & o_num_bond = atom->num_bond[o_local_index];
      int & h_num_bond = atom->num_bond[h_local_index];
      int & p_num_angle = atom->num_angle[p_local_index];
      int & o_num_angle = atom->num_angle[o_local_index];
      int & h_num_angle = atom->num_angle[h_local_index];

      // The oxygen section for the new bond
      atom->bond_type[o_local_index][o_num_bond] = typeOH;
      atom->bond_atom[o_local_index][o_num_bond] = h_global_index;
      o_num_bond++;

      // Just one bond is added for each added hydrogen atom
      nbonds_change++;

      // The oxygen section for the new angle
      atom->angle_type[o_local_index][o_num_angle] = typePOH;
      atom->angle_atom1[o_local_index][o_num_angle] = p_global_index;
      atom->angle_atom2[o_local_index][o_num_angle] = o_global_index;
      atom->angle_atom3[o_local_index][o_num_angle] = h_global_index;
      o_num_angle++;

      // Just one angle is added for each extra hydrogen atom
      nangles_change++;
      
      if (force->newton_bond) continue; // If the newton flag is on just half of the topology information is required
      
      // The hydrogen section for the new bond
      atom->bond_type[h_local_index][0] = typeOH;
      atom->bond_atom[h_local_index][0] = o_global_index;
      h_num_bond++;

      // The phosphate section for the new angle
      atom->angle_type[p_local_index][o_num_angle] = typePOH;
      atom->angle_atom1[p_local_index][o_num_angle] = p_global_index;
      atom->angle_atom2[p_local_index][o_num_angle] = o_global_index;
      atom->angle_atom3[p_local_index][o_num_angle] = h_global_index;
      p_num_angle++;

      // The Oxygen section for the new angle
      atom->angle_type[o_local_index][1] = typePOH;
      atom->angle_atom1[o_local_index][i] = p_global_index;
      atom->angle_atom2[o_local_index][i] = o_global_index;
      atom->angle_atom3[o_local_index][i] = h_global_index;
      o_num_angle++;

      // The Hydrogen section for the new angle
      atom->angle_type[h_local_index][i] = typePOH;
      atom->angle_atom1[h_local_index][i] = p_global_index;
      atom->angle_atom2[h_local_index][i] = o_global_index;
      atom->angle_atom3[h_local_index][i] = h_global_index;
      h_num_angle = 1;
   }


   // Summing up the total change in the natoms, nbonds and nangles
   MPI_Allreduce(&natoms_change,&natoms_change_total,1,MPI_INT,MPI_SUM,world);
   MPI_Allreduce(&nbonds_change,&nbonds_change_total,1,MPI_INT,MPI_SUM,world);
   MPI_Allreduce(&nangles_change,&nangles_change_total,1,MPI_INT,MPI_SUM,world);

   // Updating the number of atoms, bonds and angles
   atom->natoms += natoms_change_total;
   atom->nbonds += nbonds_change_total;
   atom->nangles += nangles_change_total;
}

/* --------------------------------------------------------------------------
   The function to remove specific number of hydrogen atoms with the local
   ids in the hAtoms[3] array
   
   
   HERE I AM NOT SURE ABOUT THE NEWTON FLAG, SHOULD THE CODE BE DIFFERENT 
   WHEN FORCE->NEWTON_BOND????
   -------------------------------------------------------------------------- */

void FixAdaptiveProtonation::remove_hydrogens(const int i, int _numHs_2_del, const int oAtoms[3], const int hAtoms[3])
{
   int nlocal = atom->nlocal;
   int numHs_2_del = _numHs_2_del;


   int p_local_index = i;
   int p_global_index = atom->tag[p_local_index];
   int & p_num_bond = atom->num_bond[p_local_index]; // it is here just for completeness, we are not adding any new bonds to the P atom here.
   int & p_num_angle = atom->num_angle[p_local_index];

   // Updating the topology for the P atom removing the angle involving the Hydrogen atoms
   for (int k = 0; k < p_num_angle; k++)
   {
      if (atom->angle_atom1[p_local_index][k] == atom->tag[hAtoms[0]] ||
          atom->angle_atom3[p_local_index][k] == atom->tag[hAtoms[0]] ||
          atom->angle_atom1[p_local_index][k] == atom->tag[hAtoms[1]] ||
          atom->angle_atom3[p_local_index][k] == atom->tag[hAtoms[1]] ||
          atom->angle_atom1[p_local_index][k] == atom->tag[hAtoms[2]] ||
          atom->angle_atom3[p_local_index][k] == atom->tag[hAtoms[2]]) {
         atom->angle_type[p_local_index][k] = atom->angle_type[p_local_index][p_num_angle-1];
         atom->angle_atom1[p_local_index][k] = atom->angle_atom1[p_local_index][p_num_angle-1];
         atom->angle_atom2[p_local_index][k] = atom->angle_atom2[p_local_index][p_num_angle-1];
         atom->angle_atom3[p_local_index][k] = atom->angle_atom3[p_local_index][p_num_angle-1];
         p_num_angle--;
         nangles_change--;
      }
   }

   

   // Updating the topology for the O atom removing the bond involving the deleted atoms
   for (int k = 0; k < 3; k++) {
      int o_local_index = oAtoms[k];
      int o_global_index = atom->tag[o_local_index];
      int o_num_bond = atom->num_bond[o_local_index];

      for (int l = 0; l < o_num_bond; l++) {
         if (atom->bond_atom[o_local_index][l] == atom->tag[hAtoms[0]] ||
             atom->bond_atom[o_local_index][l] == atom->tag[hAtoms[1]] ||
             atom->bond_atom[o_local_index][l] == atom->tag[hAtoms[2]]) {
               atom->bond_type[o_local_index][l] = atom->bond_type[o_local_index][o_num_bond-1];
               atom->bond_atom[o_local_index][l] = atom->bond_type[o_local_index][o_num_bond-1];
               o_num_bond--;
               nbonds_change--;
             }
      }
   }
   
   int k = 0;
   while(k < nlocal && numHs_2_del) {
      if (hAtoms[0] == k ||
          hAtoms[1] == k ||
          hAtoms[2] == k) {
         atom->avec->copy(nlocal-1, k , 1);
         nlocal--;
         numHs_2_del--;
         natoms_change--;
      }
      k++;
   }

   if (numHs_2_del)
      error->warning(FLERR,"There is not enough hydrogen atoms for the fix adaptive protonation to delete, this should not have happen!!!");
      
   // Summing up the change in natoms, nbonds and nangles from different procs
   MPI_Allreduce(&natoms_change,&natoms_change_total,1,MPI_INT,MPI_SUM,world);
   MPI_Allreduce(&nbonds_change,&nbonds_change_total,1,MPI_INT,MPI_SUM,world);
   MPI_Allreduce(&nangles_change,&nangles_change_total,1,MPI_INT,MPI_SUM,world);
   
   // Updating the total number of atoms, bonds and angles
   atom->natoms += natoms_change_total;
   atom->nbonds += nbonds_change_total;
   atom->nangles += nangles_change_total;  
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

