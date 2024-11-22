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
ComputeStyle(compute_GFF_constant_pH,ComputeGFFConstantPH);
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
  void init() override;
  void setup() override;
  void compute_vector() override;
  void compute_peratom() override {}; // I just wanted LAMMPS to consider this as peratom compute so the peratom energies be tallied in this timestep.
  int pack_reverse_comm(int , int , double *) override;
  void unpack_reverse_comm(int , int *, double *) override;
  
  
  
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


  double compute_Hs();
  void modify_q();
  void update_lmp();
  void compute_q_total();
  double compute_epair();
};

}    // namespace LAMMPS_NS

#endif
#endif
