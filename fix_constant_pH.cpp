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
/* ---v0.01.12----- */

#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif

#include "fix.h"
#include "fix_constant_pH.h"

#include "atom.h"
#include "atom_masks.h"
#include "pair.h"
#include "error.h"

#include "force.h"
#include "group.h"
#include "memory.h"
#include "timer.h"
#include "comm.h"
#include "kspace.h"
#include "update.h"
#include "math_const.h"
#include "modify.h"

#include <cstring>
#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixConstantPH::FixConstantPH(LAMMPS *lmp, int narg, char **arg):
  Fix(lmp, narg, arg)
{
  if (narg < 9) utils::missing_cmd_args(FLERR,"fix constant_pH", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Illegal fix constant_pH every value {}", nevery);
  typeH = utils::inumeric(FLERR,arg[4],false,lmp);
  if (typeH > atom->ntypes) error->all(FLERR,"Illegal fix constant_pH atom type {}",typeH); 
  typeHW = utils::inumeric(FLERR,arg[5],false,lmp);
  if (typeHW > atom->ntypes) error->all(FLERR,"Illegal fix constant_pH atom type {}",typeHW);
  // For hydronium the initial charges are qO=-0.833, qH1=0.611, qH2=0.611, qH3=0.611 (based on TIP3P water model)

	
  pK = utils::numeric(FLERR, arg[6], false, lmp);
  pH = utils::numeric(FLERR, arg[7], false, lmp);
  T = utils::numeric(FLERR, arg[8], false, lmp);
  
  pstyle = utils::strdup(arg[9]);
  


  qHs = 0.0;
  qHWs = 0.612;

  GFF_flag = false;
  print_Udwp_flag = false;
  int iarg = 10;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "GFF") == 0)
    {
       GFF_flag = true;
       fp = fopen(arg[iarg+1],"r");
       if (fp == nullptr)
         error->one(FLERR, "Cannot find fix constant_pH the GFF correction file {}",arg[iarg+1]);
       iarg = iarg + 2;
    }
    else if ((strcmp(arg[iarg],"Qs") == 0))
    {
       qHs = utils::numeric(FLERR, arg[iarg+1],false,lmp);
       qHWs = utils::numeric(FLERR, arg[iarg+2],false,lmp);
       iarg = iarg + 3;
    }
    else if ((strcmp(arg[iarg],"Print_Udwp") == 0))
    {
	print_Udwp_flag = true;
	Udwp_fp = fopen(arg[iarg+1],"w");
	if (Udwp_fp == nullptr) 
	    error->one(FLERR, "Cannot find fix constant_pH the Udwp debugging file {}",arg[iarg+1]);
	iarg = iarg + 2;
    }
    else
       error->all(FLERR, "Unknown fix constant_pH keyword: {}", arg[iarg]);
  }
  
  fixgpu = nullptr;

}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   #ifdef DEBUG
       std::cout << "Releasing epsilon" << std::endl;
   #endif
   memory->destroy(epsilon_init);
   #ifdef DEBUG
       std::cout << "Releasing GFF" << std::endl;
   #endif
   if (GFF) memory->destroy(GFF);

   #ifdef DEBUG
       std::cout << "Releasing pparam1" << std::endl;
   #endif
   delete [] pparam1;
   #ifdef DEBUG
       std::cout << "Releasing pstyle" << std::endl;
   #endif
   delete [] pstyle;

   if (fp && (comm->me == 0)) fclose(fp);
   if (Udwp_fp && (comm->me == 0)) fclose(Udwp_fp); // We should never reach that point as this file is writting just at the setup stage and then it will be closed
}

/* ----------------------------------------------------------------------

   ---------------------------------------------------------------------- */

int FixConstantPH::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE; // the 1st half of the velocity verlet algorithm --> update half step v_lambda and update lambda
  mask |= POST_FORCE; // calculate the lambda 
  mask |= FINAL_INTEGRATE; // the 2nd half of the velocity verlet algorithm --> update t
  return mask;	
}

/* ----------------------------------------------------------------------
   Setup
   ---------------------------------------------------------------------- */
   
void FixConstantPH::init()
{
   std::map<std::string, std::string> pair_params;
   
   pair_params["lj/cut/soft/omp"] = "lambda";
   pair_params["lj/cut/coul/cut/soft/gpu"] = "lambda";
   pair_params["lj/cut/coul/cut/soft/omp"] = "lambda";
   pair_params["lj/cut/coul/long/soft"] = "lambda";
   pair_params["lj/cut/coul/long/soft/gpu"] = "lambda";
   pair_params["lj/cut/coul/long/soft/omp"] = "lambda";
   pair_params["lj/cut/tip4p/long/soft"] = "lambda";
   pair_params["lj/cut/tip4p/long/soft/omp"] = "lambda";
   pair_params["lj/charmm/coul/long/soft"] = "lambda";
   pair_params["lj/charmm/coul/long/soft/omp"] = "lambda";
   pair_params["lj/class2/soft"] = "lambda";
   pair_params["lj/class2/coul/cut/soft"] = "lambda";
   pair_params["lj/class2/coul/long/soft"] = "lambda";
   pair_params["coul/cut/soft"] = "lambda";
   pair_params["coul/cut/soft/omp"] = "lambda";
   pair_params["coul/long/soft"] = "lambda";
   pair_params["coul/long/soft/omp"] = "lambda";
   pair_params["tip4p/long/soft"] = "lambda";
   pair_params["tip4p/long/soft/omp"] = "lambda";
   pair_params["morse/soft"] = "lambda";
   
   pair_params["lj/charmm/coul/charmm"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/gpu"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/intel"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/kk"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/omp"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/implicit"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/implicit/kk"] = "lj14_1";
   pair_params["lj/charmm/coul/charmm/implicit/omp"] = "lj14_1";
   pair_params["lj/charmm/coul/long"] = "lj14_1";
   pair_params["lj/charmm/coul/long/gpu"] = "lj14_1";
   pair_params["lj/charmm/coul/long/intel"] = "lj14_1";
   pair_params["lj/charmm/coul/long/kk"] = "lj14_1";
   pair_params["lj/charmm/coul/long/opt"] = "lj14_1";
   pair_params["lj/charmm/coul/long/omp"] = "lj14_1";
   pair_params["lj/charmm/coul/msm"] = "lj14_1";
   pair_params["lj/charmm/coul/msm/omp"] = "lj14_1";
   pair_params["lj/charmmfsw/coul/charmmfsh"] = "lj14_1";
   pair_params["lj/charmmfsw/coul/long"] = "lj14_1";
   pair_params["lj/charmmfsw/coul/long/kk"] = "lj14_1";

   if (pair_params.find(pstyle) == pair_params.end())
      error->all(FLERR,"The pair style {} is not currently supported in fix constant_pH",pstyle);
   
   pparam1 = new char[pair_params[pstyle].length()+1];
   std::strcpy(pparam1,pair_params[pstyle].c_str());
   
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::setup(int /*vflag*/)
{
   // default values from Donnini, Ullmann, J Chem Theory Comput 2016 - Table S2
    w = 50.0;
    s = 0.3;
    h = 10.0;
    k = 6.267;
    a = 0.05130;
    b = 0.001411;
    r = 21.428;
    m = 0.1078;
    d = 5.0;
    // m_lambda = 20u taken from https://www.mpinat.mpg.de/627830/usage
    m_lambda = 200;


    pair1 = nullptr;
  
    if (lmp->suffix_enable)
        pair1 = force->pair_match(std::string(pstyle)+"/"+lmp->suffix,1);
    if (pair1 == nullptr)
        pair1 = force->pair_match(pstyle,1); // I need to define the pstyle variable
    void *ptr1 = pair1->extract(pparam1,pdim1);
    if (ptr1 == nullptr)
        error->all(FLERR,"Fix constant_pH pair style {} was not found",pstyle);
    if (pdim1 != 2)
         error->all(FLERR,"Pair style parameter {} is not compatible with fix constant_pH", pparam1);
         
    epsilon = (double **) ptr1;
    
    int ntypes = atom->ntypes;
    memory->create(epsilon_init,ntypes+1,ntypes+1,"constant_pH:epsilon_init");

    // I am not sure about the limits of these two loops, please double check them
    for (int i = 0; i < ntypes+1; i++)
        for (int j = i; j < ntypes+1; j++)
             epsilon_init[i][j] = epsilon[i][j];


	
    int * type = atom->type;
    int nlocal = atom->nlocal;
    int * nums = new int[2];
    int * nums_local = new int[2];
    nums_local[0] = 0;
    nums_local[1] = 0;
    for (int i = 0; i < nlocal; i++)
    {
        if (type[i] == typeH)
	    nums_local[0]++;
	else if (type[i] == typeHW)
	    nums_local[1]++;
    }

    MPI_Allreduce(nums_local,nums,2,MPI_INT,MPI_SUM,world);
    num_Hs = nums[0];
    num_HWs = nums[1];

    if (num_HWs != 3*num_Hs) 
	error->warning(FLERR,"In the fix constant_pH the number of hydrogen atoms of the hydronium group should be three times the number of titrable groups");
    delete [] nums_local;
    delete [] nums;

	
    GFF_lambda = 0.0;
    if (GFF_flag)
	init_GFF();

    if (print_Udwp_flag)
	print_Udwp();
	
    fixgpu = modify->get_fix_by_id("package_gpu");
}

/* ----------------------------------------------------------------------
   1st half of the velocity verlet algorithm for the lambda
   ---------------------------------------------------------------------- */

void FixConstantPH::initial_integrate(int /*vflag*/)
{
   if (update->ntimestep % nevery == 0)
   {
      update_v_lambda();
      update_lambda();
   }
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::post_force(int vflag)
{
   if (update->ntimestep % nevery == 0) {
      compute_Hs<0>();
      calculate_df();
      calculate_dU();
      update_a_lambda();
      compute_q_total();
   }
   /* The force on hydrogens must be updated at every step otherwise at 
      steps at this fix is not active the pH would be very low and there
      will be a jump in pH in nevery steps                               */
   compute_Hs<1>();
}

/* ----------------------------------------------------------------------
   The 2nd half of the velocity verlet algorithm for the lambda parameter
   ---------------------------------------------------------------------- */

void FixConstantPH::post_integrate()
{
   if (update->ntimestep % nevery == 0)
       update_v_lambda();
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::calculate_df()
{
   f = 1.0/(1+exp(-50*(lambda-0.5)));
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
   U4 = 0.5*w*(1-erff(r*(lambda+m)));
   U5 = 0.5*w*(1+erff(r*(lambda-1-m)));
   dU1 = -((lambda-1-b)/(a*a))*U1;
   dU2 = -((lambda+b)/(a*a))*U2;
   dU3 = -((lambda-0.5)/(s*s))*U3;
   dU4 = -0.5*w*r*2*exp(-r*r*(lambda+m)*(lambda+m))/sqrt(M_PI);
   dU5 = 0.5*w*r*2*exp(-r*r*(lambda-1-m)*(lambda-1-m))/sqrt(M_PI);

    U =  U1 +  U2 +  U3 +  U4 +  U5;
   dU = dU1 + dU2 + dU3 + dU4 + dU5;
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::print_Udwp()
{
    double lambda_backup = this->lambda;
    int n_points = 100;
	
    lambda = -0.5;
    double dlambda = 2.0/(double)n_points;

    fprintf(Udwp_fp,"Lambda,U,dU\n");
    for (int i = 0; i <= n_points; i++) {
	calculate_dU();
        fprintf(Udwp_fp,"%f,%f,%f\n",lambda,U,dU);
	lambda += dlambda;
    }
    fclose(Udwp_fp);
    lambda = lambda_backup;
}

/* ----------------------------------------------------------------------- */

template <int stage>
void FixConstantPH::compute_Hs()
{
   if (stage == 0)
   {
      allocate_storage();
      backup_restore_qfev<1>();      // backup charge, force, energy, virial array values
      modify_epsilon_q(0.0); //should define a change_parameters(const int);
      update_lmp(); // update the lammps force and virial values
      HA = compute_epair(); // I need to define my own version using compute pe/atom // Check if HA is for lambda=0 
      backup_restore_qfev<-1>();        // restore charge, force, energy, virial array values
      modify_epsilon_q(1.0); //should define a change_parameters(const double);
      update_lmp();
      HB = compute_epair(); 
      backup_restore_qfev<-1>();      // restore charge, force, energy, virial array values
      deallocate_storage();
   }
   if (stage == 1)
   {
      modify_epsilon_q(lambda); //should define a change_parameters(const double);
      update_lmp();
   }
}

/* -------------------------------------------------ntypesntypes---------------------
   manage storage for charge, force, energy, virial arrays
   taken from src/FEP/compute_fep.cpp
------------------------------------------------------------------------- */

void FixConstantPH::allocate_storage()
{
  /* It should be nmax since in the case that 
     the newton flag is on the force in the 
     ghost atoms also must be update and the 
     nmax contains the maximum number of nlocal 
     and nghost atoms.
  */
  int nlocal = atom->nmax; 
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
   Forward-reverse copy function to be used in backup_restore_qfev()
   ---------------------------------------------------------------------- */

template  <int direction>
void FixConstantPH::forward_reverse_copy(double &a,double &b)
{
   if (direction == 1) a = b;
   if (direction == -1) b = a;
}

template  <int direction>
void FixConstantPH::forward_reverse_copy(double* a,double* b, int i)
{
   if (direction == 1) a[i] = b[i];
   if (direction == -1) b[i] = a[i];
}

template  <int direction>
void FixConstantPH::forward_reverse_copy(double** a,double** b, int i, int j)
{
   if (direction == 1) a[i][j] = b[i][j];
   if (direction == -1) b[i][j] = a[i][j];
}

/* ----------------------------------------------------------------------
   backup and restore arrays with charge, force, energy, virial
   taken from src/FEP/compute_fep.cpp
   backup ==> direction == 1
   restore ==> direction == -1
------------------------------------------------------------------------- */

template <int direction>
void FixConstantPH::backup_restore_qfev()
{
  int i;


  int nall = atom->nlocal + atom->nghost;
  int natom = atom->nlocal;
  if (force->newton || force->kspace->tip4pflag) natom += atom->nghost;

  double **f = atom->f;
  for (i = 0; i < natom; i++)
    for (int j = 0 ; j < 3; j++)
       forward_reverse_copy<direction>(f_orig,f,i,j);
  

  double *q = atom->q;
  for (int i = 0; i < natom; i++)
     forward_reverse_copy<direction>(q_orig,q,i);

  forward_reverse_copy<direction>(eng_vdwl_orig,force->pair->eng_vdwl);
  forward_reverse_copy<direction>(eng_coul_orig,force->pair->eng_coul);

  for (int i = 0; i < 6; i++)
	  forward_reverse_copy<direction>(pvirial_orig ,force->pair->virial, i);

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++) forward_reverse_copy<direction>(peatom_orig,peatom, i);
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++)
      for (int j = 0; j < 6; j++)
      forward_reverse_copy<direction>(pvatom_orig,pvatom,i,j);
  }

  if (force->kspace) {
     forward_reverse_copy<direction>(energy_orig,force->kspace->energy);
     for (int j = 0; j < 6; j++)
         forward_reverse_copy<direction>(kvirial_orig,force->kspace->virial,j);
	  
     if (update->eflag_atom) {
        double *keatom = force->kspace->eatom;
        for (i = 0; i < natom; i++) forward_reverse_copy<direction>(keatom_orig,keatom,i);
     }
     if (update->vflag_atom) {
        double **kvatom = force->kspace->vatom;
        for (i = 0; i < natom; i++) 
	  for (int j = 0; j < 6; j++)
             forward_reverse_copy<direction>(kvatom_orig,kvatom,i,j);
     }
  }
}

/* --------------------------------------------------------------

   -------------------------------------------------------------- */
   
void FixConstantPH::modify_epsilon_q(const double& scale)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;
    double * q = atom->q;

    // not sure about the range of these two loops
    /*for (int i = 0; i < ntypes + 1; i++)
	for (int j = i; j < ntypes + 1; j++)
	    if (type[i] == typeH || type[j] == typeH)
	    	epsilon[i][j] = epsilon_init[i][j] * scale;
    */

    for (int i = 0; i < nlocal; i++)
    {
        if (type[i] == typeH)
	    q[i] = qHs + scale;
	if (type[i] == typeHW)
	    q[i] = qHWs + (-scale) * (double) num_Hs/ (double) num_HWs;	
     }
	
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
   while(fgets(line,sizeof(line), fp) != NULL && i < GFF_size)
      {
          double _lambda, _GFF;
          line[strcspn(line,"\n")] = 0;
	  char * token = strtok(line, ",");
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
      GFF_lambda = GFF[i-1][1] + ((GFF[i][1]-GFF[i-1][1])/(GFF[i][0]-GFF[i-1][0]))*(lambda - GFF[i-1][0]);
   if (i == GFF_size - 1)
   {
      error->warning(FLERR,"Warning lambda of {} in Fix constant_pH out of the range, it usually should not happen",lambda);
      GFF_lambda = GFF[i][1] + ((GFF[i][1]-GFF[i-1][1])/(GFF[i][0]-GFF[i-1][0]))*(lambda - GFF[i][0]);
   }
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::update_a_lambda()
{
   if (GFF_flag) calculate_GFF();
   double NA = 6.022*1e23;
   double kj2kcal = 0.239006;
   double kT = force->boltz * T;
   df = 1.0;
   f = 1.0;
   double  f_lambda = -(HB-HA + kj2kcal*df*kT*log(10)*(pK-pH) + kj2kcal*dU - GFF_lambda);
   this->a_lambda = f_lambda / m_lambda;
   #ifdef DEBUG
	std::cout << "The a_lambda and f_lambda are :" << a_lambda << "," << f_lambda << std::endl;
   #endif

   double  H_lambda = (1-lambda)*HA + lambda*HB + kj2kcal*f*kT*log(10)*(pK-pH) + kj2kcal*U + (m_lambda/2.0)*(v_lambda*v_lambda); // This might not be needed. May be I need to tally this into energies.
   // I might need to use the leap-frog integrator and so this function might need to be in other functions than postforce()
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::update_v_lambda()
{
   double dt_lambda = update->dt;
   this->v_lambda += 0.5*this->a_lambda*dt_lambda;
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::update_lambda()
{
   double dt_lambda = update->dt;
   this->lambda += this->v_lambda * dt_lambda;
}

   
/* --------------------------------------------------------------------- */

void FixConstantPH::compute_q_total()
{
   double * q = atom->q;
   double nlocal = atom->nlocal;
   double q_local = 0.0;
   double q_total;
   double tolerance = 0.001;

   for (int i = 0; i <nlocal; i++)
       q_local += q[i];

    MPI_Allreduce(&q_local,&q_total,1,MPI_DOUBLE,MPI_SUM,world);

    if ((q_total >= tolerance || q_total <= -tolerance) && comm->me == 0)
    	error->warning(FLERR,"q_total in fix constant-pH is non-zero: {}",q_total);
}

/* --------------------------------------------------------------------- */

double FixConstantPH::compute_epair()
{
   //if (update->eflag_global != update->ntimestep)
   //   error->all(FLERR,"Energy was not tallied on the needed timestep");

   int natoms = atom->natoms;

   double energy_local = 0.0;
   double energy = 0.0;
   if (force->pair) energy_local += (force->pair->eng_vdwl + force->pair->eng_coul);

   /* As the bond, angle, dihedral and improper energies 
      do not change with the espilon, we do not need to 
      include them in the energy. We are interested in 
      their difference afterall */

   MPI_Allreduce(&energy_local,&energy,1,MPI_DOUBLE,MPI_SUM,world);
   energy /= static_cast<double> (natoms); // To convert to kcal/mol the total energy must be devided by the number of atoms
   return energy;
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
