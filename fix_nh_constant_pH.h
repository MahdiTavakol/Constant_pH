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

namespace LAMMPS_NS {

class FixNHConstantPH : public FixNH {
 public:
  FixNHConstantPH(class LAMMPS *, int, char **);
  ~FixNHConstantPH() override;
  void init() override;
  void setup(int) override;
  double memory_usage() override;

 protected:
  double lambda, v_lambda, t_lambda;

  void nve_x() override;
  void nve_v() override;
  void nh_v_temp() override;
  virtual int size_restart_global();

};

}    // namespace LAMMPS_NS

#endif
