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
#include "pair.h"
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
  if (narg < 9) utils::mising_cmd_args(FLERR,"fix constant_pH", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Ilegal fix constant pH every value {}", nevery);
  igroupH = group->find(arg[4]);
  if (igroupH == -1) error->all(FLERR,"Cannot find the hydrogens group for fix constant_pH);
  groupHbit = group->bitmask(igroupH);
  igroupW = group->find(arg[5]);
  if (igroupW == -1) error->all(FLERR,"Cannot find the water group for fix constant_pH");
  if (group->count(igroupW) != 3) 
     error->all(FLERR, "Number of atoms in the water molecule for the fix constant_pH is {} instead of three",group->count(igroupW));
  groupWbit = group->bitmask(igroupW);
  pK = utils::numeric(FLERR, arg[6], false, lmp);
  pH = utils::numeric(FLERR, arg[7], false, lmp);
  T = utils::numeric(FLERR, arg[8], false, lmp);
   
  int iarg = 9;
  while (iarg < narg) {
    // for now, no other keywords
  }

}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::post_force(int vflag)
{
   if (update->ntimestep % nevery == 0) {
      compute_Hs();
      calculate_df();
      calculate_dU();
      integrate_lambda();
   }
   /* The force on hydrogens must be updated at every step otherwise at 
      steps at this fix is not active the pH would be very low and there
      will be a jump in pH in nevery steps                               */
   set_force(); // I am not sure with change_parameter it should be here or not!!!!
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::init()
{
   // default values from Donnini, Ullmann, J Chem Theory Comput 2016 - Table S2
    w = 200.0;
    s = 0.3;
    h = 4.0;
    k = 2.533;
    a = 0.034041;
    b = 0.005238;
    r = 16.458;
    m = 0.1507;
    d = 2.0;
    // m_lambda = 20u taken from https://www.mpinat.mpg.de/627830/usage
    m_lambda = 20;

    // It might needs to be in the setup. I am not sure
    nmax = 1;
    if (atom->nmax > nmax) {
	memory->destroy(H_atom);
    	nmax = atom->nmax;
    	memory->create(H_atom, nmax, "constant_pH:H_atom");
    }
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

void FixConstantPH::compute_Hs()
{
   double **f = atom->f;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;

  // This part should be in the init function 
   Pair * pair = nullptr;
  
   if (lmp->suffix_enable)
       pair = force->pair_match(std::string(ad->pstyle)+"/"+lmp->suffix,1);
   if (pair == nullptr)
       pair = force->pair_match(pstyle,1); // I need to define the pstyle variable
   void *ptr1 = pair->extract(pparam1,pdim1);
   if (ptr1 == nullptr)
       error->all(FLERR,"Fix adapt/fep pair style param not supported");
   if (pdim != 2)
        error->all(FLERR,"Pair style parameter {} is not compatible with fix constant_pH", pparam1);
   
   if (pdim == 2) epsilons = (double **) ptr;

   for (int i = 0; i <= ??; i++)
        for (int j = i; j <= ??; j++)
            epsilons_org[i][j] = epsilons[i][j];

  // This part should be in the init function 





  eflag = 1;
  vflag = 0;

  invoked_vector = update->ntimestep;

  if (atom->nmax > nmax) {    // reallocate working arrays if necessary
     ???? 
    ?????
  }

  backup_qfev();      // backup charge, force, energy, virial array values

  change_parameters(0.0); //should define a change_parameters(const int);

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(eflag, vflag);
    timer->stamp(Timer::PAIR);
  }
  if (chgflag && force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(eflag, vflag);
    timer->stamp(Timer::KSPACE);
  }

  // accumulate force/energy/virial from /gpu pair styles
  if (fixgpu) fixgpu->post_force(vflag);

  HA = compute_epair(); // I need to define my own version using compute pe/atom // Check if HA is for lambda=0

  change_parameters(1.0); //should define a change_parameters(const double);

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(eflag, vflag);
    timer->stamp(Timer::PAIR);
  }
  if (chgflag && force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(eflag, vflag);
    timer->stamp(Timer::KSPACE);
  }

  // accumulate force/energy/virial from /gpu pair styles
  // this is required as to empty the answer queue,
  // otherwise the force compute on the GPU in the next step would be incorrect
  if (fixgpu) fixgpu->post_force(vflag);

  HB = compute_epair();

  restore_qfev();      // restore charge, force, energy, virial array values
  restore_params();    // restore pair parameters

  change_parameters(lambda); //should define a change_parameters(const double);



}

/* --------------------------------------------------------------------------
   -------------------------------------------------------------------------- */

void FixConstantPH::modify_params()
{
    for (int i = 0; i < ??; i++)
	for (int j = i; j < ??; j++)
	    epsilons[i][j] *= lambda;
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array --> Needs to be updated at the end
------------------------------------------------------------------------- */

double ComputePEAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
