// clang-format off
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

#include "fix.h"

#include "atom.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "math_const.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixConstantPH::FixConstantPH(LAMMPS *lmp, int narg, char **arg):
  Fix(lmp, narg, arg)
{
  if (narg < ??) utils::mising_cmd_args(FLERR,"fix constant pH", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Ilegal fix constant pH every value {}", nevery);

  int iarg = 4;
  while (iarg < narg) {
    
  }


   groupHbit = group->bitmask[igroupH];
}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::post_force(int vflag)
{
   if (update->ntimestep % nevery) return;

   calculate_df();
   calculate_dU();
   integrate_lambda();
   set_force();
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::init()
{
   // set the initial lambda and v_lambda values
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::integrate_lambda()
{
   double f_lambda = -(HB-HA + df*R*T*ln(10)*(pK-pH) + dU);
   double a_lambda = f_lambda / m_lambda;
   double t_lambda = nevery*update->dt;
   double H_lambda = (1-lambda)*HA + lambda*HB + f*R*T*ln(10)*(pK-pH) + U + (m_lambda/2.0)*(v_lambda**2);
   lambda = (1.0/2.0)*(a_lambda)*(t_lambda)**2 + v_lambda*t_lambda + lambda;
   v_lambda = a_lambda * t_lambda + v_lambda;
}
/* ---------------------------------------------------------------------- */

void FixConstantPH::calculate_df()
{
   f = 1.0/(1+exp(-50*(lambda-0.5));
   df = 50*exp(-50*(lambda-0.5))/(f**2);
}

/* ----------------------------------------------------------------------- */

void FixConstantPH::calculate_dU()
{
   double U1, U2, U3, U4, U5;
   double dU1, dU2, dU3, dU4, dU5;
   U1 = -k*exp(-(lambda-1-b)*(lambda-1-b)/(2*a*a));
   U2 = -k*exp(-(lambda+b)*(lambda+b)/(2*a*a));
   U3 = d*exp(-(lambda-0.5)*(lambda-0.5)/(2*s*s));
   U4 = 0.5*w*(1-erff(r*(lambda+m));
   U5 = 0.5*w*(1+erff(r*(lambda-1-m));
   dU1 = -((lambda-1-b)/(2*a*a))*U1;
   dU2 = -((lambda+b)/(2*a*a))*U2;
   dU3 = -((lambda-0.5)/(s*s))*U3;
   dU4 = -0.5*w*r*2*exp(-r*r*(lambda+0.5)*(lambda+0.5))/sqrt(MY_PI);
   dU5 = 0.5*w*r*2*exp(-r*r*(lambda-1-m)*(lambda-1-m))/sqrt(MY_PI);

    U =  U1 +  U2 +  U3 +  U4 +  U5;
   dU = dU1 + dU2 + dU3 + dU4 + dU5;
}

/* ----------------------------------------------------------------------- */

void FixConstantPH::set_force()
{
   double **f = atom->f;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   
   for (int i = 0; i < nlocal; i++)
   {
      if (mask[i] & groupHbit)
      {
         f[i][0] *= lambda;
         f[i][1] *= lambda;
         f[i][2] *= lambda;
      }
   }
}
