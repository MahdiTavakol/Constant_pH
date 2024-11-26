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

#ifndef LMP_FIX_NH_CONSTANT_PH_H
#define LMP_FIX_NH_CONSTANT_PH_H

#include "fix.h"    // IWYU pragma: export
#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHConstantPH : public FixNH {
 public:
  FixNHConstantPH(class LAMMPS *, int, char **);
  ~FixNHConstantPH() override;
  double memory_usage() override;

 protected:

  void nve_x() override;
  void nve_v() override;
  void nh_v_temp() override;

  // lambda variables from the fix constant pH  
  Fix *fix_constant_pH;
  char *fix_constant_pH_id;
  double v_lambda, m_lambda, a_lambda, T_lambda;

};

}    // namespace LAMMPS_NS

#endif
