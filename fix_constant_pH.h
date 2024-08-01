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
FixStyle(constant_pH,FixConstantPH);
// clang-format on
#else

#ifndef LMP_FIX_CONSTANTPH_H
#define LMP_FIX_CONSTANTPH_H

#include "fix.h"


namespace LAMMPS_NS {

	class FixConstantPH: public Fix {
	public:
		FixConstantPH(class LAMMPS*, int, char**);
		~FixConstantPH() override;
		int setmask() override;
		void init() override;
		void setup(int) override;
		void post_force(int) override;
		double compute_scalar() override;
		double compute_vector(int) override;
		double memory_usage() override;
		void init_list(int, class NeighList*) override;

	private:
		// Input variables for constant values
		int igroupH, igroupW;
		int groupHbit, groupWbit;
		double pK, pH, T;
		double a, b, s, m, w, r, d;
		double m_lambda;
		double HA, HB;
		int nmax;
                double *H_atom;
		char* pparam1; // It should be epsilon or something like that.
                char* pparam2;
                int pdim1;
                int pdim2;
                double** epsilons;
                double** epsilons_org;

		void integrate_lambda();
		void compute_Hs();
		void calculate_df();
		void calculate_dU();
		void set_force();
		void modify_water();
	};

}    // namespace LAMMPS_NS


#endif
#endif
