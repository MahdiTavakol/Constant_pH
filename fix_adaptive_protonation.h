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
/* ---v0.0.00----- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(AdaptiveProtonation,FixAdaptiveProtonation);
// clang-format on
#else

#ifndef LMP_FIX_ADAPTIVEPROTONATION_H
#define LMP_FIX_ADAPTIVEPROTONATION_H

#include "fix.h"

namespace LAMMPS_NS {

	class FixMetaAR : public Fix {
	public:
		FixMetaAR(class LAMMPS*, int, char**);
		~FixMetaAR() override;
		int setmask() override;
		void init() override;
		void setup(int) override;
		void pre_exchange(int) override;
		double compute_scalar() override;
		double compute_vector(int) override;
		double memory_usage() override;
		void init_list(int, class NeighList*) override;

	private:
		// Input variables for constant values
    		double typeP, typeO, typeH, typeOw, typeHw, typeOHbond, threshold;
    

		// Input variables for variable values
		char* typePstr, *typeOstr, *typeSstr;

		//
		int typePvar, typeOvar, typeHvar;
    		int typePstyle, typeOstyle, typeSstyle;


		// Neighborlist is required for accessing neighbors
		class NeighList* list;


		// Number of neighbors
		double * numNeighbors;

		// number of added atoms
		double numaddedatoms;

		// mark 0  -> not type P
		// mark 1  -> in water 
		// mark -1 -> in solid
                // mark -2 -> on surface
		int * mark;


		void mark_protonation_deprotonation();

    		int nmax;

	};

}    // namespace LAMMPS_NS

#endif
#endif
