/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ---v0.10.00----- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(AdaptiveProtonation,FixAdaptiveProtonation);
// clang-format on
#else

#ifndef LMP_FIX_ADAPTIVE_PROTONATION_H
#define LMP_FIX_ADAPTIVE_PROTONATION_H

#include "fix.h"
#include "atom_vec.h"

namespace LAMMPS_NS {

   class FixAdaptiveProtonation : public Fix {
   public:
      FixAdaptiveProtonation(class LAMMPS*, int, char**);
      ~FixAdaptiveProtonation() override;
      int setmask() override;
      void init() override;
      void pre_exchange() override;
      double compute_scalar() override;
      double compute_vector(int) override;
      double memory_usage() override;
      void init_list(int, class NeighList*) override;


      // Getting the number of protonable molids;
      void get_n_protonable(int & n_lambdas) const
      {
         n_lambdas = this->n_protonable;
      }
      // Getting the protonable molids (The molids array must be allocated otherwise an error is produced)
      void get_protonable_molids(int * molids) const;

   protected:

      /// -----> This part is similar to the part in the fix_constant_pH.cpp, so should be a separate class
      // The input files 
      FILE* pHstructureFile3, * pHstructureFile4;

   
      // The information on the protonable species
      int pHnStructures3, pHnStructures4;
      int pHnTypes1, pHnTypes2;
      double **pH1qs, **pH2qs;
      int * typePerProtMol;
      int *protonable;

      void read_pH_structure_files();

      // <------ This part is similar to the part in the fix_constant_pH.cpp, so should be a separate class


      // I need to detect water molecules
      int typeOW;
      // The threshold for the number of neighboring water molecules
      double threshold;

      // pKa and pH values
      int pKa, pH;
      
      // How many hydrogen atoms should there be
      int req_numHs;
    

      // Input variables for variable values
      char* typePstr, *typeOstr, *typeHstr, *typeOHstr, *typePOHstr;

      //
      int typePvar, typeOvar, typeHvar, typeOHvar, typePOHvar;
      int typePstyle, typeOstyle, typeHstyle, typeOHstyle, typePOHstyle;


      // Neighborlist is required for accessing neighbors
      class NeighList* list;


      /*
       * -1 ---> NEITHER
       *  0 ---> SOLID
       *  1 ---> SOLVENT
       */
       
      int * mark;
      int * mark_local;
      int * mark_prev; // For the previous step
      int * mark_per_mol; // If one atom have mark == 1 all the atoms of that molecule should have mark == 1

      int * molecule_size; // used to average the mark for each molecule
      int * molecule_size_local;
      
      // track changes in the natoms, nbonds and nangles
      int natoms_change, natoms_change_total;
      int nbonds_change, nbonds_change_total;
      int nangles_change, nangles_change_total;



      // array to access the molids of the protonable molecules
      int * protonable_molids;
      int n_protonable;

      int nmax;
      int nmolecules;


      // Setting the same molecule id for atoms connected to each other
      void set_molecule_id();
      // Allocates the protonable_molids array based on the number of molecule_ids obtained from set_molecule_id()
      void set_protonable_molids();
      // Mark phosphate atoms for protonation/deprotonation
      void mark_protonation_deprotonation();
      // Fills out the array of the molids of the protonable molecules
      void reset_protonable_molids();
      // Add/Remove hydrogens for all the phosphates in the system based on the output of the mark_protonation_deprotonation() function
      void modify_protonable_hydrogens();
      // Add hydrogen atoms to the phosphate with a P atom having an index of i
      void add_hydrogens(const int& i, const int& req_numHs, const int oAtoms[3], const int hAtoms[3]);
      // Remove hydrogen atoms from the phosphate with a P atom having an index of i
      void remove_hydrogens(const int i, int numHs_2_del, const int oAtoms[3], const int hAtoms[3]);

   };

}    // namespace LAMMPS_NS

#endif
#endif
