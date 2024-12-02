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

/* ----------------------------------------------------------------------
   Constant pH support added by: Mahdi Tavakol (Oxford)
   v0.03.06
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NH_CONSTANT_PH_H
#define LMP_FIX_NH_CONSTANT_PH_H

#include <ctime>

#include "fix.h"    // IWYU pragma: export
#include "fix_constant_pH.h"
#include "fix_nh.h"


namespace LAMMPS_NS {

class FixNHConstantPH : public FixNH {
 public:
  FixNHConstantPH(class LAMMPS *, int, char **);
  ~FixNHConstantPH() override;
  void init() override;
  double memory_usage() override;

 protected:

  void nve_x() override;
  void nve_v() override;
  void nh_v_temp() override;
  double random_normal(double mean, double stddev);

  // lambda variables from the fix constant pH  
  FixConstantPH *fix_constant_pH;
  char *fix_constant_pH_id;
  double* x_lambdas, *v_lambdas, *a_lambdas, *m_lambdas;
  double T_lambda;
  int n_lambdas;

};

}    // namespace LAMMPS_NS

#endif
