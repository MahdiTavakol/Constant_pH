/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_temp_constant_pH.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempConstantPH::ComputeTempConstantPH(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg), fix_constant_pH_id(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute temp command");
  fix_constant_pH_id = utils::strdup(arg[3]);

  scalar_flag = vector_flag = 1;
  size_vector = 7; // I need to double check to see if the fix_nh.cpp can access the seventh element or if this element causes any problem for the fix_nh.cpp
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeTempConstantPH::~ComputeTempConstantPH()
{
  if (!copymode) delete[] vector;
  if (fix_constant_pH_id) delete[] fix_constant_pH_id;
}

/* ---------------------------------------------------------------------- */

void ComputeTempConstantPH::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();

  fix_constant_pH_id = modify->get_fix_by_id(fix_constant_pH_id);
}


/* ---------------------------------------------------------------------- */

void ComputeTempConstantPH::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp + 1; // The +1 is for the lambda dof
  dof -= extra_dof + fix_dof;
  if (dof > 0.0)
    tfactor = force->mvv2e / (dof * force->boltz);
  else
    tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempConstantPH::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }

  double a_lambda;
  fix_constant_pH->return_params(&m_lambda,&v_lambda,&a_lambda); // The return_parameters section should be implemented in the fix_constant_pH.cpp code
  
  t += v_lambda*v_lambda * m_lambda; // I need to define a method to access the v_lambda and m_lambda from the fix_constant_pH.h

  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}
