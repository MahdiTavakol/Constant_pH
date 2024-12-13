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
/* ---------- v0.05.04----------------- */
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
   Fix(lmp, narg, arg), typePstr(nullptr), typeOstr(nullptr), typeHstr(nullptr)
{
   if (narg < 11) utils::missing_cmd_args(FLERR, "fix AdaptiveProtonation", error);

   dynamic_group_allow = 0;
   scalar_flag = 1; 
   vector_flag = 1;
   size_vector = 3; // What is this?
   global_freq = 1; // What is this?
   extscalar = 1; // What is this?
   extvector = 1; // What is this?
   energy_global_flag = 1;
   virial_global_flag = 1;
   virial_peratom_flag = 1;
   respa_level_support = 1;
   ilevel_respa = 0;

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
      typeHstyle = CONSTANT;
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
   
   nmax = 1;
   memory->create(mark,nmax,"AdaptiveProtontation:mark");
}

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::~FixAdaptiveProtonation()
{
   delete[] typePstr;
   delete[] typeOstr;
   delete[] typeHstr;
   delete[] typeOHstr;

   memory->destory(mark);
}

/* --------------------------------------------------------------------------------------- */

int FixAdaptiveProtonation::setmark()
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
}

/* ---------------------------------------------------------------------------------------
    It is need to access the neighbor list
   --------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::init_list(int /*id*/, NeighList* ptr) 
{
   list = ptr;
}

/* --------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::pre_exchange(int vflag)
{
   if(update->ntimestep != next_reneighbor) return;

   if (atom->nmax > nmax)
   {
      memory->destroy(mark);
      nmax = atom->nmax;
      memory->create(mark,nmax,"AdaptiveProtontation:mark");
   }

   mark_protonation_deprotonation();
   change_protonation();  
}

/* ----------------------------------------------------------------------------------------
   getting the number of water molecules near phosphates and 
   flag them for protonation/deprotonation
   ---------------------------------------------------------------------------------------- */

void FixAdaptiveProtonation::mark_protonation_deprotonation()
{
   int *ilist, *jlist, *numneigh, **firstneight;
   int inum;
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
         continue
      }
      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; j++) {
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
   int * molecule = atom->molecule;
   int * num_bond = atom->num_bond;
   int ** bond_atom = atom->bond_atom;
   int * tag = atom->tag; // atom-id
   int * molecule_id = new int[natom];

   for (int i = 0; i < natom; i++) {
      molecule_id[i] = tag[i];
   }

   for (int i = 0; i < nlocal; i++) {
      for (int k = 0; k < num_bond[i]; k++) {
         int jtag = bond_atom[i][k]; // the tag (atom-id) of kth bonds of atom i
         int j = atom->map(jtag);
         if (j == -1) {
            error->warning(FLERR,"Bond atom missing in fix AdaptiveProtonation");
            continue;
         }
         molecule_id[i] = MIN(molecule_id[i],molecule_id[j]); // I am not sure about header for the MIN
         molecule_id[j] = molecule_id[i];
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

   for (int i = 0; i < nlocal; i++)
      molecule[i] = molecule_id[i];

   delete [] molecule_id;
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


   if (atom->avec->bonds_allow && (style == BOND || style == MULTI || style == ATOM)) {
      
      int pAtom;     // local atom id for the P atom
      int oAtoms[3]; // local atom id for the O atoms
      int hAtoms[3]; // local atom id for the H atoms

      for (int i = 0; i < nlocal; i++) {
         if (mark[i] == 0) continue;

         for (int j = 0; j < nlocal; j++)
            hlist[j] = 0;

         pAtom = i;
         int numHs = 0; 
         for (int m = 0; m < atom->num_bond[i]; m++) {
            /* 
               bond_atoms[i][m] contains the global id (tag) of the connected atom
               so there is a need for mapping the global id to the local one
            */
            int atom2 = atom->map(atom->bond_atom[i][m]);
            if (m < 3) oAtoms[m] = atom2;
            else error->warning(FLERR,"Number of atoms connected to P is {}",num_bond[i]);
            if (atom2 == -1) error->warning(FLERR,"Bond atom missing in fix AdaptiveProtonation");
            for (int n = 0; n < atom->num_bond[atom2]; n++) {
               int atom3 = atom->map(atom->bond_atom[atom2][n]);
               if (type[atom3] == typeH) 
               {
                   if (numHs < 3)
                      hAtoms[numHs++] = atom3;
                   else
                      error->warning(FLERR,"Number of protons are {}",++numHs);
                  hlist[atom3] = 1;
               }
            }

            /* 
               Considering different values of mark[i] and setting the number of hydrogen atoms
               to an appropriate value
            */
            
            if (mark[i] == -1 && numHs == 0) continue; // PO4 inside the HAp ---> Nothing to do here!
            if (mark[i] == -1 && numHs > 0) remove_hydrogens(numHs,hlist); //HPO4, H2PO4 or H3PO4 inside ---> The hydrogens must be removed
            if (mark[i] == 1 && numHs == req_numHs) continue; // HxPO4 inside the surface with the exact number of required Hydrogen atoms ---> Nothing to do here!
            if (mark[i] == 1 && numHs > req_numHs) remove_hydrogens(numHs-req_numHs,hlist); // HxPO4 in the surface or in solution with higher number of Hydrogens
            if (mark[i] == 1 && numHs < req_numHs) {
               /* HxPO4 in the surface or in solution with lower number of Hydrogens
                  I would prefer to remove all the hydrogens of this phosphate and then add
                  the required number to eliminate the possibility of having an Oxygen atom with
                  two hydrogen atoms
               */
               remove_hydrogens(numHs, hAtoms);
               add_hydrogens(i,req_numHs,oAtoms,hAtoms);
            }
         }
      }
   }
}

void FixAdaptiveProtonation::add_hydrogens(const int& i, const int& req_numHs, const int* const oAtoms, const int* const hAtoms) 
{
   double *x = atom->x;

   // adding hydrogens to this phosphate
   double coors[3][3]; // Maximum three phosphate ions
   int tags[3];
   int js[3];
               
   coors[0][0] = x[oAtoms[0]][0] + 0.9;
   coors[0][1] = x[oAtoms[0]][1];
   coors[0][2] = x[oAtoms[0]][2];
   coors[1][0] = x[oAtoms[1]][0] - 0.8;
   coors[1][1] = x[oAtoms[1]][1];
   coors[1][2] = x[oAtoms[1]][2];
   coors[2][0] = x[oAtoms[2]][0];
   coors[2][1] = x[oAtoms[2]][1] + 0.64;
   coors[2][2] = x[oAtoms[2]][2] - 0.64;


   error->warning(FLERR,"atom->natoms before adding atoms = {}",atom->natoms);
   atom->avec->create_atom(typeH,coor1);
   atom->avec->create_atom(typeH,coor2);
   atom->avec->create_atom(typeH,coor3);
   error->warning(FLERR,"atom->natoms after adding atoms = {}",atom->natoms);
               
   numaddedatoms+=2;

   tags[0] = (tag[nlocal-3] > 0)? tag[nlocal-3]: atom->natoms-3;
   tags[1] = (tag[nlocal-2] > 0)? tag[nlocal-2]: atom->natoms-2;
   tags[2] = (tag[nlocal-1] > 0)? tag[nlocal-1]: atom->natoms-1;
   atom->natom += 3; // I am not sure if this is required or not.


   // atom->avec->create_atoms itself updates the atom->nlocal value
   js[0] = atom->nlocal - 3;
   js[1] = atom->nlocal - 2;
   js[2] = atom->nlocal - 1;
               
   // Let's connect these two new Hs to the P
   if (atom->num_bond[i] == atom->bond_per_atom)
      error->one(FLERR,"Num bonds exceeded bonds per atom in fix AdaptiveProtonation");

   /*
     I guess that the size of bond_type and bond_atom is atom->bond_per_atom so there is enough space there
   */

   for (int k = 0; k < req_numHs; k++) {
      // The kth OH bond
      int o_local_index = oAtoms[k];
      int h_local_index = hAtoms[k];
      int o_global_index = atom->tag[oAtoms[k]];
      int & o_num_bond = atom->num_bond[o_local_index];
      int & h_num_bond = atom->num_bond[h_local_index];

      // The oxygen section
      atom->bond_type[o_local_index][o_num_bond] = bondOHtype;
      atom->bond_atom[o_local_index][o_num_bond] = tags[k];
      o_num_bond++;

      // The hydrogen section
      atom->bond_type[h_local_index][0] = bondOHtype;
      atom->bond_atom[h_local_index][0] = o_global_index;
      h_num_bond++;
   }
}

void FixAdaptiveProtonation::remove_hydrogens(int& numHs_2_del, int hAtoms[3])
{
   int nlocal = atom->nlocal;
   int k = 0;
   
   while(k < nlocal && numHs_2_del) {
      if (hAtoms[0] == k ||
          hAtoms[1] == k ||
          hAtoms[2] == k) {
         avec->copy(nlocal-1, k , 1);
         nlocal--;
         numHs_2_del--;
      }
      k++;
   }

   if (numHs_2_del)
      error->warning(FLERR,"There is not enough hydrogen atoms for the fix adaptive protonation to delete, it should not have happen!");
}
