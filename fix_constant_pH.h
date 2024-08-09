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

#ifdef FIX_CLASS
// clang-format off
FixStyle(constant_pH,FixConstantPH);
// clang-format on
#else

#ifndef LMP_FIX_CONSTANTPH_H
#define LMP_FIX_CONSTANTPH_H

#include "fix.h"
#include "pair.h"


namespace LAMMPS_NS {

  class FixConstantPH: public Fix {
     public:
	FixConstantPH(class LAMMPS*, int, char**);
	~FixConstantPH() override;
	int setmask() override;
	void init() override;
	void initial_integrate(int) override;
	void post_force(int) override;
	void post_integrate() override;
	double memory_usage() override;

     private:
	// Input variables for constant values
	int typeH, typeHW;
	double pK, pH, T;

	double a, b, s, m, w, r, d, k, h;
	double HA, HB;
	double U, dU;
	
	// Pair style parameters
	char * pstyle, * pparam1;
	Pair * pair1;
	int pdim1;

	// Lambda dynamics
	double lambda, v_lambda, a_lambda, m_lambda;	

	// Protonation and hydronium group parameters
	double qHs, qHWs;
	int num_Hs, num_HWs;

	// The smoothing function 
	double f, df;

	// Parameters for the forcefield modifiction term
        bool GFF_flag;
	FILE *fp;
	double **GFF;
	int GFF_size;
	double GFF_lambda;


	class Fix *fixgpu;

	// This is just a pointer to the non-bonded interaction parameters and does not have any allocated memory
	// This should not be deallocated since the original pointer will be deallocated later on by the LAMMPS
	double ** epsilon;
	// _init is the initial value of hydrogen atoms properties which is multiplied by lambda at each step
	double **epsilon_init;

        // _org is for value of parameters before the update_lmp() with modified parameters act on them
  	double *q_orig;
 	double **f_orig;
  	double eng_vdwl_orig, eng_coul_orig;
  	double pvirial_orig[6];
  	double *peatom_orig, **pvatom_orig;
 	double energy_orig;
 	double kvirial_orig[6];
	double *keatom_orig, **kvatom_orig;


        template<int stage>
	void compute_Hs();
	void calculate_df();
	void calculate_dU();
	void integrate_lambda();
	void allocate_storage();
	void deallocate_storage();
        template < int direction>
	void forward_reverse_copy(double& a, double& b);
        template < int direction>
	void forward_reverse_copy(double* a, double* b, int i);
        template < int direction>
	void forward_reverse_copy(double** a, double** b, int i, int j);
	template <int direction>
	void backup_restore_qfev();
	void init_GFF();
	void calculate_GFF();
	void modify_epsilon_q(const double& scale);
	void modify_water();
	void update_lmp();
	double compute_epair();
	void update_a_lambda();
	void update_v_lambda();
	void update_lambda();
	};

}    // namespace LAMMPS_NS


#endif
#endif
