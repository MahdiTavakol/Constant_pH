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
/* ---------- v0.05.03----------------- */
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
   if (narg < 9) utils::missing_cmd_args(FLERR, "fix AdaptiveProtonation", error);

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
      error->all(FLERR,"The variable input for typeOH bond in AdaptiveProtonation is not supported yet!");
   }
   else
   {
      typeOHbond = utils::numeric(FLERR,arg[7],false,lmp);
   }
   if (utils::strmatch(arg[8], "^v_")) {
      error->all(FLERR,"The variable input for threshold in AdaptiveProtonation is not supported yet!");
   }
   else
   {
      threshold = utils::numeric(FLERR,arg[8],false,lmp);
   }


   // optional args



   // set up reneighboring from src/fix_evaporate.cpp
   force_reneighbor = 1;
   next_reneighbor = (update->ntimestep/nevery) * nevery + nevery;

   // the number of added atoms
   numaddedatoms = 0;

   nmax = 1;
   memory->create(mark,nmax,"AdaptiveProtontation:mark");
}

/* --------------------------------------------------------------------------------------- */

FixAdaptiveProtonation::~FixAdaptiveProtonation()
{
   delete[] typePstr;
   delete[] typeOstr;
   delete[] typeHstr;

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
         int j;
         int jtag = bond_atom[i][k]; // the tag (atom-id) of kth bonds of atom i
         for (int j = 0; j < natom; j++) {
            if (tag[j] == jtag)
               break;
         }
         if (j == natom) {
            error->warning(FLERR,"fix adaptive protonation cannot find the bonded atom to atom with id of {}",tag[i]);
            continue;
         }
         molecule_id[j] = jtag;
         // The question is that what is what the molecule id of jtag atom is 
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

void FixAdaptiveProtonation::change_protonation()
{
   int nlocal = atom->nlocal;
   tagint *tag = atom->tag;
   int * type = atom->type;
   double **x = atom->x;
   AtomVec *avec = atom->avec;
   int * hlist = new int(nlocal);

   for (int i = 0; i < nlocal; i++)
      hlist[i] = 0;

   if (atom->avec->bonds_allow && 
      (style == BOND || style == MULTI || style == ATOM)) {
      int * num_bond = atom->num_bond;

      for (int i = 0; i < nlocal; i++) {
         if (mark[i] == 0) continue;
         int pAtom, * oAtoms, * hAtoms; 
         oAtoms = new int[3];
         hAtoms = new int[3];
         pAtom = i;
         int numHs = 0; 
         for (int m = 0; m < num_bond[i]; m++) {
            int atom2 = atom->map(atom->bond_atom[i][m]);
            if (m < 3) oAtoms[i] = atom2;
            else error->warning(FLERR,"Number of atoms connected to P is {}",num_bond[i]);
            if (atom2 == -1) error->warning(FLERR,"Bond atom missing in fix AdaptiveProtonation");
            for (int n = 0; n < num_bond[atom2]; n++) {
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
            if (mark[i] == 1  && numHs == 3) continue; // H3PO4 in solution. Nothing to do here!  
            if (mark[i] == -1 && numHs == 0) continue; // PO4 out of water. Nothing to do here!
            // This has definitly came from the HAp surface
            if (mark[i] == 1 && numHs == 0) {
               // adding hydrogens to this phosphate
               // For the moment we just add two hydrogens
               double * coor1, * coor2;
               coor1[0] = x[i][0] + 0.9;
               coor1[1] = x[i][1];
               coor1[2] = x[i][2];
               coor2[0] = x[i][0] - 0.8;
               coor2[1] = x[i][1];
               coor2[2] = x[i][2];
               atom->avec->create_atom(typeH,coor1);
               atom->avec->create_atom(typeH,coor2);
               int j1 = atom->natoms + 1;
               int j2 = atom->natoms + 2;
               atom->natom += 1;
               // Let's connect these two new Hs to the P
               if (num_bond[i] == atom->bond_per_atom)
                  error->one(FLERR,"Num bonds exceeded bonds per atom in fix AdaptiveProtonation");
               bond_type[i][num_bond[i]] = bondOHtype;
               bond_atom[i][num_bond[i]] = tag[j1];
               num_bond[i]++;
               if (num_bond[i] == atom->bond_per_atom)
                  error->one(FLERR,"Num bonds exceeded bonds per atom in fix AdaptiveProtonation");
               bond_type[i][num_bond[i]] = bondOHtype;
               bond_atom[i][num_bond[i]] = tag[j2];
               num_bond[i]++;
            }
            // This has come from the water
            if (mark[i] == -1 && numHs > 0) {
               while (numHs) {
                  while (k < nlocal) {
                     if (hlist[k] == 1) {
                       avec->copy(nlocal - 1, k, 1);
                       nlocal--;
                       numHs--;
                     } else
                    k++;
                  } 
               }
            }
         }
         delete [] oAtoms;
         delete [] hAtoms;
      }
   }

   delete [] hlist;
}
