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
/* ---v0.00.3----- */

#include "fix.h"

#include "atom.h"
#include "atom_masks.h"
#include "pair.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "timer.h"
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
  if (nevery < 0) error->all(FLERR,"Illegal fix constant_pH every value {}", nevery);
  typeH = utils::inumeric(FLERR,arg[4],flase,lmp);
  if (typeH > atom->ntypes) error->all(FLERR,"Illegal fix constant_pH atom type {}",typeH); 
  igroupW = group->find(arg[5]);
  if (igroupW == -1) error->all(FLERR,"Cannot find the hydronium ions hydrogen groups for fix constant_pH");
  // For hydronium the initial charges are qO=-0.834, qH1=0.611, qH2=0.612, qH3=0.611 (based on TIP3P water model)
  // if m is the number of igroupW atoms, m-1 atoms in the igroupW should be given a charge of lambda/3.0 
  // and then the last atom must be given the charge of m*lambda/3 - (m-1)*lambda/3. This is done to count for rounding errors
  groupWbit = group->bitmask(igroupW);
	
  pK = utils::numeric(FLERR, arg[6], false, lmp);
  pH = utils::numeric(FLERR, arg[7], false, lmp);
  T = utils::numeric(FLERR, arg[8], false, lmp);

  GFF_flag = false;
  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "GFF") == 0)
    {
       GFF_flag = true;
       fp = fopen(arg[iarg+1],"r");
       if (fp == nullptr)
         error->one(FLERR, "Cannot find fix constant_pH the GFF correction file {}",arg[iarg+1]);
       iarg = iarg + 2;
    }
    if ((strcmp(arg[iarg],"qs") == 0)
    {
       qHW = utils::numeric(FLERR, arg[iarg+1]);
       iarg = iarg + 1;
    }
    else
       error->all(FLERR, "Unknown fix constant_pH keyword: {}", arg[iarg]);
  }

}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   memory->destroy(epsilon_init); 
   memory->destroy(GFF);

   if (fp && (comm->me == 0)) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::post_force(int vflag)
{
   if (update->ntimestep % nevery == 0) {
      compute_Hs<0>();
      calculate_df();
      calculate_dU();
      integrate_lambda();
   }
   /* The force on hydrogens must be updated at every step otherwise at 
      steps at this fix is not active the pH would be very low and there
      will be a jump in pH in nevery steps                               */
   compute_Hs<1>();
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


    pair1 = nullptr;
  
    if (lmp->suffix_enable)
        pair1 = force->pair_match(std::string(pstyle)+"/"+lmp->suffix,1);
    if (pair1 == nullptr)
        pair1 = force->pair_match(pstyle,1); // I need to define the pstyle variable
    void *ptr1 = pair1->extract(pparam1,pdim1);
    if (ptr1 == nullptr)
        error->all(FLERR,"Fix constant_pH pair style param not supported");
    if (pdim1 != 2)
         error->all(FLERR,"Pair style parameter {} is not compatible with fix constant_pH", pparam1);
   
    if (pdim1 == 2) double ** epsilon = (double **) ptr1;

    int ntypes = atom->ntypes;
    memory->create(epsilon_init,ntypes+1,ntypes+1,"constant_pH:epsilon_init");

    // I am not sure about the limits of these two loops, please double check them
    for (int i = 0; i <= ntypes+1; i++)
         for (int j = i; j <= ntypes+1; j++)
             epsilon_init[i][j] = epsilon[i][j];

	
    int * type = atom->type;
    int nlocal = atom->nlocal;
    int num_Hs_local = 0;
    for (int i = 0; i < nlocal; i++)
	if (type[i] == typeH)
	    num_Hs_local++;

    MPI_Allreduce(&num_Hs_local,&num_Hs,MPI_INT,MPI_SUM,world);
    num_HWs = group->count(groupWbit);

    if (num_HWs != 3*num_Hs) 
	error->warning(FLERR,"In the fix constant_pH the number of hydrogen atoms of the hydronium group should be three times the number of titrable groups");
	
    GFF_lambda = 0.0;
    if (GFF_flag)
	init_GFF();
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::integrate_lambda()
{
   if (GFF_flag) calculate_GFF();
   double  f_lambda = -(HB-HA + df*R*T*ln(10)*(pK-pH) + dU - GFF_lambda);
   double  a_lambda = f_lambda / m_lambda;
   double dt_lambda = update->dt;
   double  H_lambda = (1-lambda)*HA + lambda*HB + f*R*T*ln(10)*(pK-pH) + U + (m_lambda/2.0)*(v_lambda**2); // This might not be needed. May be I need to tally this into energies.
   // I might need to use the leap-frog integrator and so this function might need to be in other functions than postforce()
   lambda += (1.0/2.0)*(a_lambda)*(dt_lambda)*(dt_lambda) + v_lambda*dt_lambda;
   v_lambda += a_lambda * dt_lambda;
}
/* ---------------------------------------------------------------------- */

void FixConstantPH::calculate_df()
{
   f = 1.0/(1+exp(-50*(lambda-0.5));
   df = 50*exp(-50*(lambda-0.5))/(f*f);
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

template <int stage>
void FixConstantPH::compute_Hs()
{
   if (stage == 0)
   {
      if (atom->nmax > nmax) {    // reallocate working arrays if necessary
         deallocate_storage();
         allocate_storage();
      }
      backup_qfev();      // backup charge, force, energy, virial array values
      modify_params(0.0); //should define a change_parameters(const int);
      update_lmp(); // update the lammps force and virial values
      HA = compute_epair(); // I need to define my own version using compute pe/atom // Check if HA is for lambda=0
      restore_qfev();      // restore charge, force, energy, virial array values
      restore_params();    // restore pair parameters and charge values
      modify_params(1.0); //should define a change_parameters(const double);
      update_lmp();
      HB = compute_epair();
      restore_qfev();      // restore charge, force, energy, virial array values
      restore_params();    // restore pair parameters and charge values
   }
   if (stage == 1)
   {
      modify_params(lambda); //should define a change_parameters(const double);
      update_lmp();
   }
}

/* ----------------------------------------------------------------------
   manage storage for charge, force, energy, virial arrays
   taken from src/FEP/compute_fep.cpp
------------------------------------------------------------------------- */

void FixConstantPH::allocate_storage()
{
  int nlocal = atom->nlocal; // nlocal is fine since before the next atom exchange between MPI ranks these variables will be destroyed.
  memory->create(f_orig, nlocal, 3, "constant_pH:f_orig");
  memory->create(q_orig, nlocal, "constant_pH:q_orig");
  memory->create(peatom_orig, nlocal, "constant_pH:peatom_orig");
  memory->create(pvatom_orig, nlocal, 6, "constant_pH:pvatom_orig");
  if (force->kspace) {
     memory->create(keatom_orig, nlocal, "constant_pH:keatom_orig");
     memory->create(kvatom_orig, nlocal, 6, "constant_pH:kvatom_orig");
  }
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::deallocate_storage()
{
  memory->destroy(q_orig);
  memory->destroy(f_orig);
  memory->destroy(peatom_orig);
  memory->destroy(pvatom_orig);
  memory->destroy(keatom_orig);
  memory->destroy(kvatom_orig);

  f_orig = nullptr;
  peatom_orig = keatom_orig = nullptr;
  pvatom_orig = kvatom_orig = nullptr;
}

/* ----------------------------------------------------------------------
   backup and restore arrays with charge, force, energy, virial
   taken from src/FEP/compute_fep.cpp
------------------------------------------------------------------------- */

void ComputeFEP::backup_qfev()
{
  int i;

  int nall = atom->nlocal + atom->nghost;
  int natom = atom->nlocal;
  if (force->newton || force->kspace->tip4pflag) natom += atom->nghost;

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f_orig[i][0] = f[i][0];
    f_orig[i][1] = f[i][1];
    f_orig[i][2] = f[i][2];
  }

  double *q = atom->q;
  for (int i = 0; i < natom; i++)
     q_orig[i] = q[i];

  eng_vdwl_orig = force->pair->eng_vdwl;
  eng_coul_orig = force->pair->eng_coul;

  pvirial_orig[0] = force->pair->virial[0];
  pvirial_orig[1] = force->pair->virial[1];
  pvirial_orig[2] = force->pair->virial[2];
  pvirial_orig[3] = force->pair->virial[3];
  pvirial_orig[4] = force->pair->virial[4];
  pvirial_orig[5] = force->pair->virial[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++) peatom_orig[i] = peatom[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom_orig[i][0] = pvatom[i][0];
      pvatom_orig[i][1] = pvatom[i][1];
      pvatom_orig[i][2] = pvatom[i][2];
      pvatom_orig[i][3] = pvatom[i][3];
      pvatom_orig[i][4] = pvatom[i][4];
      pvatom_orig[i][5] = pvatom[i][5];
    }
  }

  if (force->kspace) {
     energy_orig = force->kspace->energy;
     kvirial_orig[0] = force->kspace->virial[0];
     kvirial_orig[1] = force->kspace->virial[1];
     kvirial_orig[2] = force->kspace->virial[2];
     kvirial_orig[3] = force->kspace->virial[3];
     kvirial_orig[4] = force->kspace->virial[4];
     kvirial_orig[5] = force->kspace->virial[5];

     if (update->eflag_atom) {
        double *keatom = force->kspace->eatom;
        for (i = 0; i < natom; i++) keatom_orig[i] = keatom[i];
     }
     if (update->vflag_atom) {
        double **kvatom = force->kspace->vatom;
        for (i = 0; i < natom; i++) {
          kvatom_orig[i][0] = kvatom[i][0];
          kvatom_orig[i][1] = kvatom[i][1];
          kvatom_orig[i][2] = kvatom[i][2];
          kvatom_orig[i][3] = kvatom[i][3];
          kvatom_orig[i][4] = kvatom[i][4];
          kvatom_orig[i][5] = kvatom[i][5];
        }
     }
  }
}

/* --------------------------------------------------------------------------
   -------------------------------------------------------------------------- */

void FixConstantPH::modify_params(const double& scale)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;

    // not sure about the range of these two loops
    for (int i = 0; i < ntypes + 1; i++)
	for (int j = i; j < ntypes + 1; j++)
	    if (type[i] == typeH || type[j] == typeH)
	    	epsilon[i][j] = epsilon_init[i][j] * scale;

    int numWH = group->count(igroupW);
    for (int i = 0; i < nlocal; i++)
    {
        if (type[i] == typeH)
	    q[i] = scale;
	if (mask[i] & groupWbit)
	    q[i] = qHW + (-scale) * (double) numHs/ (double) numHWs;	
     }
	    
	    
    // I need to add bond, angle, dihedral, improper and charge parameters
	
}

/* ----------------------------------------------------------------------
   modify force and kspace in lammps according
   ---------------------------------------------------------------------- */

void FixConstantPH::update_lmp() {
   int eflag = 1;
   int vflag = 0;
   timer->stamp();
   if (force->pair && force->pair->compute_flag) {
     force->pair->compute(eflag, vflag);
     timer->stamp(Timer::PAIR);
   }
   if (force->kspace && force->kspace->compute_flag) {
     force->kspace->compute(eflag, vflag);
     timer->stamp(Timer::KSPACE);
   }

   // accumulate force/energy/virial from /gpu pair styles
   if (fixgpu) fixgpu->post_force(vflag);
}

/* ---------------------------------------------------------------------
   Read the data file containing the term deltaGFF in equation 2 of 
   https://pubs.acs.org/doi/full/10.1021/acs.jctc.5b01160
   --------------------------------------------------------------------- */

void FixConstantPH::init_GFF()
{
   char line[100];
   fgets(line,sizeof(line),fp);
   line[strcspn(line,"\n")] = 0;
   GFF_size = atoi(line);
   memory->create(GFF,GFF_size,2,"constant_pH:GFF");
   int i = 0;
   while(fgets(line,sizeof(line), fp) != NULL && i < GFF_array_size)
      {
          line[strcspn(line,"\n")] = 0;
	  token = strtok(line, ",");
	  if (token != NULL) 
	      _lambda = atof(token);
	  else 
	      error->one(FLERR,"The GFF correction file in the fix constant_pH has a wrong format!");
	  token = strtok(line, ",");
	  if (token != NULL)
	      _GFF = atof(token);
	  else
	      error->one(FLERR,"The GFF correction file in the fix constant_pH has a wrong format!");
	  GFF[i][0] = _lambda;
	  GFF[i][1] = _GFF;
	  i++;   
      }	
    if (fp && (comm->me == 0)) fclose(fp);
}

/* ---------------------------------------------------------------------
   Add forcefield correction term deltaGFF in equation 2 of
   https://pubs.acs.org/doi/full/10.1021/acs.jctc.5b01160
   --------------------------------------------------------------------- */

void FixConstantPH::calculate_GFF()
{
   int i = 0;
   while (i < GFF_size && GFF[i][0] < lambda) i++;

   if (i == 0)
   {
      error->warning(FLERR,"Warning lambda of {} in Fix constant_pH out of the range, it usually should not happen",lambda);
      GFF_lambda = GFF[0][1] + ((GFF[1][1]-GFF[0][1])/(GFF[1][0]-GFF[0][0]))*(lambda - GFF[0][0]);
   }
   if (i > 0 && i < GFF_size - 1)
      GFF_lambda = GFF[i-1] + ((GFF[i][1]-GFF[i-1][1])/(GFF[i][0]-GFF[i-1][0]))*(lambda - GFF[i-1][0]);
   if (i == GFF_size - 1)
   {
      error->warning(FLERR,"Warning lambda of {} in Fix constant_pH out of the range, it usually should not happen",lambda);
      GFF_lambda = GFF[i] + ((GFF[i][1]-GFF[i-1][1])/(GFF[i][0]-GFF[i-1][0]))*(lambda - GFF[i][0]);
   }
}
   

/* ----------------------------------------------------------------------
   memory usage of local atom-based array --> Needs to be updated at the end
   ---------------------------------------------------------------------- */

double FixConstantPH::memory_usage()
{
  int nmax = atom->nmax;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  double pair_bytes = sizeof(Pair);
  double GFF_bytes = 2.0 * GFF_size * sizeof(double);
  double epsilon_init_bytes = (double) (ntypes + 1) * (double) (ntypes + 1) * sizeof(double);
  double q_orig_bytes = (double) nlocal * sizeof(double);
  double f_orig_bytes = (double) nlocal * 3.0 * sizeof(double);
  double peatom_orig_bytes = (double) nlocal * sizeof(double);
  double pvatom_orig_bytes = (double) nlocal * 6.0 * sizeof(double);
  double keatom_orig_bytes = (double) nlocal * sizeof(double);
  double kvatom_orig_bytes = (double) nlocal * 6.0 * sizeof(double);
  double bytes = pair_bytes + GFF_bytes + epsilon_init_bytes + \
	         q_orig_bytes + f_orig_bytes + peatom_orig_bytes + \
                 pvatom_orig_bytes + keatom_orig_bytes + kvatom_orig_bytes;
  return bytes;
}
