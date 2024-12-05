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
	void setup(int) override;
	void initial_integrate(int) override;
        double compute_vector(int) override;
	double memory_usage() override;

        // Functions for accessing or reseting the lambda dynamics parameters
        void return_nparams(int& _n_params) const;
        void return_params(double* const _x_lambdas, double* const _v_lambdas, 
                           double* const _a_lambdas, double* const _m_lambdas) const;
        void reset_params(const double* const _x_lambdas, const double* const _v_lambdas, 
                          const double* const _a_lambdas, const double* const _m_lambdas);

     private:
	// Sturcture files
        FILE *pHStructureFile;

	// Atom types and charges that change due to protonation
        int pHnTypes;
        double *pH1qs, *pH2qs;
        int * typePerProtMol;
        int * protonable;


	// Input variables for constant values
	int typeHW;
	double pK, pH, T;

	double a, b, s, m, w, r, d, k, h;
	double HA, HB;
	double U, dU;
	
	// Pair style parameters
        // I am not sure why I do not release the pstyle
	char * pstyle, * pparam1;
	Pair * pair1;
	int pdim1;

	// Lambda dynamics
	double lambda, v_lambda, a_lambda, m_lambda, H_lambda;

        // Lambda arrays
        double * lambdas, * v_lambdas, * a_lambdas, * m_lambdas;
        int n_lambdas;

        // The protonable groups
        int *protonable_molecule_ids;
         

	// Protonation and hydronium group parameters
	double qHs, qHWs;
	int num_HWs, num_prots;

	// The smoothing function 
	double f, df;

	// Parameters for the forcefield modifiction term
        bool GFF_flag;
	FILE *fp;
	double **GFF;
	int GFF_size;
	double GFF_lambda;

        // Parameters for printing the Udwp
        bool print_Udwp_flag;
        FILE *Udwp_fp;
        void print_Udwp();
        
        // Parameters for thermostating the lambda variable
        double etha_lambda;
        double Q_lambda;
        double T_lambda;
        double a_etha_v_ratio_lambda;

        // The q_total used to calculate the HW charges
        double q_total;
        void compute_q_total();
        void check_q_total();


	class Fix *fixgpu;

	// This is just a pointer to the non-bonded interaction parameters and does not have any allocated memory
	// This should not be deallocated since the original pointer will be deallocated later on by the LAMMPS
	double **epsilon;
	// _init is the initial value of hydrogen atoms properties which is multiplied by lambda at each step
	double **epsilon_init;


        int nmax;

        // These pointers are allocated and deallocated through allocate_storage() and deallocate_storage() functions
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
        void calculate_num_prot_num_HWs();
        void read_pH_structure_files();
        void restore_epsilon();
	void calculate_dq();
	void calculate_df();
	void calculate_dU();
	void calculate_T_lambda();
	void initialize_v_lambda(const double _T_lambda);
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
	void modify_q(const double& scale);
	void update_lmp();
        void compute_f_lambda_charge_interpolation();
	double compute_epair();
	void update_a_lambda();
	};

}    // namespace LAMMPS_NS


#endif
#endif
