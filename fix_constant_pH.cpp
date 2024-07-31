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
   set_force(); 
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

/* --------------------------------------------------------------------------
     taken from src/compute_pe_atom.cpp
   -------------------------------------------------------------------------- */

void FixConstantPH::compute_Hs()
{
  int i;

  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR, "Per-atom energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(H_atom);
    nmax = atom->nmax;
    memory->create(H_atom, nmax, "constant_pH:H_atom");
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set

  int nlocal = atom->nlocal;
  int npair = nlocal;
  int nbond = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton_bond) nbond += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  if (force->kspace && force->kspace->tip4pflag) nkspace += atom->nghost;

  // clear local energy array

  for (i = 0; i < ntotal; i++) H_atom[i] = 0.0;

  // add in per-atom contributions from each force

  if (pairflag && force->pair && force->pair->compute_flag) {
    double *eatom = force->pair->eatom;
    for (i = 0; i < npair; i++) H_atom[i] += eatom[i];
  }

  if (bondflag && force->bond) {
    double *eatom = force->bond->eatom;
    for (i = 0; i < nbond; i++) H_atom[i] += eatom[i];
  }

  if (angleflag && force->angle) {
    double *eatom = force->angle->eatom;
    for (i = 0; i < nbond; i++) H_atom[i] += eatom[i];
  }

  if (dihedralflag && force->dihedral) {
    double *eatom = force->dihedral->eatom;
    for (i = 0; i < nbond; i++) H_atom[i] += eatom[i];
  }

  if (improperflag && force->improper) {
    double *eatom = force->improper->eatom;
    for (i = 0; i < nbond; i++) H_atom[i] += eatom[i];
  }

  if (kspaceflag && force->kspace && force->kspace->compute_flag) {
    double *eatom = force->kspace->eatom;
    for (i = 0; i < nkspace; i++) H_atom[i] += eatom[i];
  }

  // add in per-atom contributions from relevant fixes
  // always only for owned atoms, not ghost

  if (fixflag && modify->n_energy_atom) modify->energy_atom(nlocal, H_atom);

  // communicate ghost energy between neighbor procs

  if (force->newton || (force->kspace && force->kspace->tip4pflag)) comm->reverse_comm(this);

  // calculate HA and HB

  int *mask = atom->mask;

  double HA_local = 0.0;
  double HB_local = 0.0;
  double * Hs_local = new double[2];
  double * Hs_all = new double[2];

  for (i = 0; i < nlocal; i++) {
    HA_local += H_atom[i];
    if (!(mask[i] & groupHbit)) HB_local += H_atom[i];
  }
  // You need to consider the water molecule here 


  Hs[0] = HA_local;
  Hs[1] = HB_local;
  // The final magic!
  MPI_Allreduce(&Hs_local, &Hs_all, 2 , MPI_DOUBLE, MPI_SUM, world);

  HA = Hs_all[0];
  HB = Hs_all[1];
  delete [] Hs_local;
  delete [] Hs_all;
}

/* May be I need to set the size of the reverse communication through 
   variable comm_reverse inherited from the Fix class. I have no idea
   how does this communication occur                                  */
/* ---------------------------------------------------------------------- */

int FixConstantPH::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = H_atom[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    H_atom[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePEAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
