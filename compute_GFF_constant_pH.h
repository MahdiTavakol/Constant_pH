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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(constant_pH/GFF,ComputeGFFConstantPH);
// clang-format on
#else

#ifndef COMPUTE_GFF_CONSTANT_PH_H
#define COMPUTE_GFF_CONSTANT_PH_H

#include "compute.h"
#include "pair.h"

namespace LAMMPS_NS {

class ComputeGFFConstantPH : public Compute {
 public:
  ComputeGFFConstantPH(class LAMMPS *, int, char **);
  ~ComputeGFFConstantPH() override;
  void setup() override;
  void init() override;
  void compute_vector() override;
  void compute_peratom() override {}; // I just wanted LAMMPS to consider this as peratom compute so the peratom energies be tallied in this timestep.
  
  
  
 private:
  // Sturcture files
  FILE *pHStructureFile;

  // dlambda for the calculation of dU/dlambda
  double lambda, dlambda;



  // Atom types and charges that change due to protonation
  int pHnTypes;
  double *pH1qs, *pH2qs;
  int * typePerProtMol;
  int * protonable;

  // Charge difference between structure 1 and structure 2
  double dq;


  // The protonium hydrogen atoms
  int typeHW;
  int num_HWs;
  double qHWs;

  class Fix *fixgpu;

  // HA ==> H(lambda-dlambda), HB ==> H(lambda+dlambda), HC ==> H(lambda), dH_dLambda ==> (HB-HA)/(2*dlambda)
  double HA, HB, HC, dH_dLambda;
  
  // _org is for value of parameters before the update_lmp() with modified parameters act on them
  double *q_orig;
  double **f_orig;
  double eng_vdwl_orig, eng_coul_orig;
  double pvirial_orig[6];
  double *peatom_orig, **pvatom_orig;
  double energy_orig;
  double kvirial_orig[6];
  double *keatom_orig, **kvatom_orig;

  int nmax;
  bigint natoms;

  // Energy of each atom computed for states A and B
  double *energy_peratom;
  
  
  

  void allocate_storage();
  void deallocate_storage();

  template  <int direction>
  static void forward_reverse_copy(double &,double &);
  template  <int direction>
  static void forward_reverse_copy(double* ,double* , int );
  template  <int direction>
  static void forward_reverse_copy(double** ,double** , int , int );
  template <int direction>
  void backup_restore_qfev();
  
  
  // Reading the charges for the deprotonated and protonated states
  void read_pH_structure_files();
  // Compute the HA, HB, HC and dH_dLambda values
  void compute_Hs();
  // Modifying the lambda value
  void modify_lambda(const double& scale);
  // Recalculating the energy values
  void update_lmp();
  // Calculating the total pair energy of the system
  double compute_epair();
  // The change in the q as a result of protonation
  void calculate_dq();
  // Checking if we have enough HW atoms to neutralize the system
  void check_num_HWs();
  // Debugging the total charge of the system
  void compute_q_total();
};

}    // namespace LAMMPS_NS

#endif
#endif
