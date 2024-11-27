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
#include "fix.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempConstantPH::ComputeTempConstantPH(LAMMPS *lmp, int narg, char **arg) : ComputeTemp(lmp, narg, arg), 
fix_constant_pH_id(nullptr), x_lambdas(nullptr), v_lambdas(nullptr), a_lambdas(nullptr), m_lambdas(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute temp constant pH command");
  fix_constant_pH_id = utils::strdup(arg[3]);

  n_lambdas = 1;
  int iarg = 4;
  while (iarg < narg) {
     if (!strcmp(arg[iarg], "n_lambdas")) {
        n_lambdas = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        iarg += 2;
     } else error->all(FLERR,"Illegal compute temp constant pH command");
  }

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
  if (x_lambdas) delete[] x_lambdas;
  if (v_lambdas) delete[] v_lambdas;
  if (a_lambdas) delete[] a_lambdas;
  if (m_lambdas) delete[] m_lambdas;
}

/* ---------------------------------------------------------------------- */

void ComputeTempConstantPH::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();

  fix_constant_pH = static_cast<FixConstantPH*>(modify->get_fix_by_id(fix_constant_pH_id));

  x_lambdas = new double[n_lambdas];
  v_lambdas = new double[n_lambdas];
  a_lambdas = new double[n_lambdas];
  m_lambdas = new double[n_lambdas];
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
  int _n_lambdas;

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

  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  fix_constant_pH->return_nparams(_n_lambdas);
  if (n_lambdas != _n_lambdas)
     error->all(FLERR,"The n_lambdas parameter in the compute temperature constant pH is not the same as the n_lambdas in the fix constant pH: {},{}",n_lambdas,_n_lambdas);
  fix_constant_pH->return_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas); // The return_parameters section should be implemented in the fix_constant_pH.cpp code

  for (int i = 0; i < n_lambdas; i++)
     scalar += v_lambdas[i]*v_lambdas[i] * m_lambdas[i];

   
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}
