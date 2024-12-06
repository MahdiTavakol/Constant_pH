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
/* ---v0.04.02----- */

#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif
#include <random>

#include "fix.h"
#include "fix_constant_pH.h"

#include "atom.h"
#include "atom_masks.h"
#include "error.h"

#include "force.h"
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "timer.h"
#include "comm.h"
#include "kspace.h"
#include "update.h"
#include "math_const.h"
#include "modify.h"

#include <cstring>
#include <string>
#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixConstantPH::FixConstantPH(LAMMPS *lmp, int narg, char **arg): Fix(lmp, narg, arg), 
       lambdas(nullptr), v_lambdas(nullptr), a_lambdas(nullptr), m_lambdas(nullptr), H_lambdas(nullptr),
       protonable(nullptr), typePerProtMol(nullptr), pH1qs(nullptr), pH2qs(nullptr),
       HAs(nullptr), HBs(nullptr), GFF_lambdas(nullptr), molids(nullptr),
       fs(nullptr), dfs(nullptr), Us(nullptr), dUs(nullptr),
       GFF(nullptr),
       q_orig(nullptr), f_orig(nullptr),
       peatom_orig(nullptr), pvatom_orig(nullptr), keatom_orig(nullptr), kvatom_orig(nullptr)
{
  if (narg < 9) utils::missing_cmd_args(FLERR,"fix constant_pH", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Illegal fix constant_pH every value {}", nevery);
  // Reading the file that contains the charges before and after protonation/deprotonation
  if (comm->me == 0) {
      pHStructureFile = fopen(arg[4],"r"); // The command reads the file the type and charge of each atom before and after protonation
      if (pHStructureFile == nullptr)
         error->all(FLERR,"Unable to open the file");
  }
  // Hydronium ion hydrogen atoms
  typeHW = utils::inumeric(FLERR,arg[5],false,lmp);
  if (typeHW > atom->ntypes) error->all(FLERR,"Illegal fix constant_pH atom type {}",typeHW);
  // For hydronium the initial charges are qO=-0.833, qH1=0.611, qH2=0.611, qH3=0.611 (based on TIP3P water model)

	
  pK = utils::numeric(FLERR, arg[6], false, lmp);
  pH = utils::numeric(FLERR, arg[7], false, lmp);
  T = utils::numeric(FLERR, arg[8], false, lmp);
  


  qHs = 0.0;
  qHWs = 0.0; //0.278;

  GFF_flag = false;
  print_Udwp_flag = false;
  n_lambdas = 1;
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
    else if ((strcmp(arg[iarg],"molids") == 0)) {
	n_lambdas = utils::numeric(FLERR,arg[iarg+1],false,lmp);
	memory->create(molids,n_lambdas,"constant_pH:lambdas");
	iarg+=2;
	for (int i = 0; i < n_lambdas; i++) {
	    molids[i] = utils::numeric(FLERR,arg[iarg],false,lmp);
	    iarg++;
	}
    }
    else
       error->all(FLERR, "Unknown fix constant_pH keyword: {}", arg[iarg]);
   }
  
   fixgpu = nullptr;
  
  
  
   array_flag = 1;
   size_array_rows = 10;
   size_array_cols = n_lambdas;

}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   if (lambdas)   memory->destroy(lambdas);
   if (v_lambdas) memory->destroy(v_lambdas);
   if (a_lambdas) memory->destroy(a_lambdas);
   if (m_lambdas) memory->destroy(m_lambdas);
   if (H_lambdas) memory->destroy(H_lambdas);
   if (protonable) memory->destroy(protonable);
   if (typePerProtMol) memory->destroy(typePerProtMol);
   if (pH1qs) memory->destroy(pH1qs);
   if (pH2qs) memory->destroy(pH2qs);
   if (HAs) memory->destroy(HAs);
   if (HBs) memory->destroy(HBs);
   if (GFF_lambdas) memory->destroy(GFF_lambdas);
   if (molids) memory->destroy(molids);
   if (fs) memory->destroy(fs);
   if (dfs) memory->destroy(dfs);
   if (Us) memory->destroy(Us);
   if (dUs) memory->destroy(dUs);

   if (GFF) memory->destroy(GFF);
   

   deallocate_storage();

   if (fp && (comm->me == 0)) fclose(fp);
   if (Udwp_fp && (comm->me == 0)) fclose(Udwp_fp); // We should never reach that point as this file is writting just at the setup stage and then it will be closed
}

/* ----------------------------------------------------------------------

   ---------------------------------------------------------------------- */

int FixConstantPH::setmask()
{
   int mask = 0;
   mask |= INITIAL_INTEGRATE; // Calculates the a_lambda
   return mask;	
}

/* ----------------------------------------------------------------------
   Setup
   ---------------------------------------------------------------------- */
   
void FixConstantPH::init()
{
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::setup(int /*vflag*/)
{
   // default values from Donnini, Ullmann, J Chem Theory Comput 2016 - Table S2
    w = 200;
    s = 0.3;
    h = 4;
    k = 2.553;
    a = 0.03401;
    b = 0.005238;
    r = 16.458; 
    m = 0.1507;
    d = 2.0;


	
    // Reading the structure of protonable states before and after protonation.
    read_pH_structure_files();


    // Checking if we have enough hydronium ions to neutralize the system
    calculate_num_prot_num_HWs();
	
    fixgpu = modify->get_fix_by_id("package_gpu");





    memory->create(lambdas,n_lambdas,"constant_pH:lambdas");
    memory->create(v_lambdas,n_lambdas,"constant_pH:v_lambdas");
    memory->create(a_lambdas,n_lambdas,"constant_pH:a_lambdas");
    memory->create(m_lambdas,n_lambdas,"constant_pH:m_lambdas");
    memory->create(H_lambdas,n_lambdas,"constant_pH:H_lambdas");

    memory->create(HAs,n_lambdas,"constant_pH:HAs");
    memory->create(HBs,n_lambdas,"constant_pH:HBs");
    memory->create(GFF_lambdas,n_lambdas,"constant_pH:GFF_lambdas");
    memory->create(fs,n_lambdas,"constant_pH:fs");
    memory->create(dfs,n_lambdas,"constant_pH:df");
    memory->create(Us,n_lambdas,"constant_pH:Us");
    memory->create(dUs,n_lambdas,"constant_pH:dUs");
    
    
    for (int j = 0; j < n_lambdas; j++)
        GFF_lambdas[j] = 0.0;
    if (GFF_flag)
	init_GFF();

    if (print_Udwp_flag)
	print_Udwp();
	
    // This would not work in the initialize section as the m_lambda has not been set yet!
    initialize_v_lambda(this->T);
	

    for (int i = 0; i < n_lambdas; i++) {
	lambdas[i] = 0.0;
	v_lambdas[i] = 0.0;
	a_lambdas[i] = 0.0;
	m_lambdas[i] = 20.0; // m_lambda == 20.0u taken from https://www.mpinat.mpg.de/627830/usage
    }
	
    nmax = atom->nmax;
    allocate_storage();

}

/* ----------------------------------------------------------------------
   This part calculates the acceleration of the lambdas parameter
   which is obtained from the force acting on it
   ---------------------------------------------------------------------- */

void FixConstantPH::initial_integrate(int /*vflag*/)
{
   compute_Hs<-1>();
   calculate_dfs();
   calculate_dUs();
   update_a_lambda();
   compute_Hs<1>();
   compute_q_total();
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::update_a_lambda()
{
   if (GFF_flag) calculate_GFFs();
   double NA = 6.022*1e23;
   double kj2kcal = 0.239006;
   double kT = force->boltz * T;

   //df = 1.0;
   //f = 1.0;

   for (int i = 0; i < n_lambdas; i++) {
	double  f_lambda = -(HBs[i]-HAs[i] - dfs[i]*kT*log(10)*(pK-pH) + kj2kcal*dUs[i] - GFF_lambdas[i]); // I'm not sure about the sign of the df*kT*log(10)*(pK-pH) 
	this->a_lambdas[i] = f_lambda /m_lambdas[i]; // 4.184*0.0001*f_lambda / m_lambda;
	// I am not sure about the sign of the f*kT*log(10)*(pK-pH)
        this->H_lambdas[i] = (1-lambdas[i])*HAs[i] + lambdas[i]*HBs[i] - fs[i]*kT*log(10*(pK-pH)) + kj2kcal*Us[i] + (m_lambdas[i]/2.0)*(v_lambdas[i]*v_lambdas[i]); // This might not be needed. May be I need to tally this into energies.
        // I might need to use the leap-frog integrator and so this function might need to be in other functions than postforce()
   }	
}
	
/* ----------------------------------------------------------------------- */

template <int stage>
void FixConstantPH::compute_Hs()
{
   if (stage == -1)
   {
      if (nmax < atom->nmax)
      {
	  nmax = atom->nmax;
          allocate_storage();
	  deallocate_storage();
      }
      for (int j = 0; j < n_lambdas; j++) {
	   double* lambdas_j = new double[n_lambdas];
	   std::fill(lambdas_j,lambdas_j+n_lambdas,0.0);
	   backup_restore_qfev<1>();
	   lambdas_j[j] = 0.0;
	   modify_qs(lambdas_j);
	   update_lmp();
	   HAs[j] = compute_epair();
           backup_restore_qfev<-1>();
	   lambdas_j[j] = 1.0;
	   modify_qs(lambdas_j);
	   HBs[j] = compute_epair();
	   backup_restore_qfev<-1>();
           delete [] lambdas_j;
      }
   }
   if (stage == 1)
   {
      modify_qs(lambdas); //should define a change_parameters(const double);
      //update_lmp(); This update_lmp() might not work here since I am not sure about the correct values for the eflag and vflag variables... Anyway, the epsilon and charge values have been updated according to the pH value and lammps will do the rest
   }
}

/* ----------------------------------------------------------------------
   returns the number of the lambda parameters
  ----------------------------------------------------------------------- */

void FixConstantPH::return_nparams(int& _n_params) const
{
    _n_params = this->n_lambdas;
}

/* ----------------------------------------------------------------------
   returns the x_lambdas, v_lambdas, ....
   The memories for these should be allocated before hands
   ---------------------------------------------------------------------- */

void FixConstantPH::return_params(double* const _x_lambdas, double* const _v_lambdas, 
                           double* const _a_lambdas, double* const _m_lambdas) const 
{
    for (int i = 0; i < n_lambdas; i++) {
	_x_lambdas[i] = lambdas[i];
	_v_lambdas[i] = v_lambdas[i];
	_a_lambdas[i] = a_lambdas[i];
	_m_lambdas[i] = m_lambdas[i];
    }
}

/* ---------------------------------------------------------------------
    This one just returns the value of T_lambda
   --------------------------------------------------------------------- */
   
void FixConstantPH::return_T_lambda(double& _T_lambda)
{
    calculate_T_lambda();
    _T_lambda = this->T_lambda;
}

/* ---------------------------------------------------------------------
    sets the values of the x_lambdas, v_lambdas, ... possibly by the intergrating
    fix styles
    --------------------------------------------------------------------- */

void FixConstantPH::reset_params(const double* const _x_lambdas, const double* const _v_lambdas, 
                          const double* const _a_lambdas, const double* const _m_lambdas)
{
    for (int i = 0; i < n_lambdas; i++) {
	lambdas[i] = _x_lambdas[i];
	v_lambdas[i] = _v_lambdas[i];
	a_lambdas[i] = _a_lambdas[i];
	m_lambdas[i] = _m_lambdas[i];
    }
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::read_pH_structure_files()
{
   /* File format
    * Comment 
    * pHnTypes
    * type1,  number of type1 atoms in the protonable molecule, qBeforeProtonation, qAfterProtonation
    * ...   U1 = -k*exp(-(lambda-1-b)*(lambda-1-b)/(2*a*a));
   U2 = -k*exp(-(lambda+b)*(lambda+b)/(2*a*a));
   U3 = d*exp(-(lambda-0.5)*(lambda-0.5)/(2*s*s));
   U4 = 0.5*w*(1-erff(r*(lambda+m)));
   U5 = 0.5*w*(1+erff(r*(lambda-1-m)));
    * ...
    */

   /*Allocating the required memory*/
   int ntypes = atom->ntypes;
   memory->create(protonable,ntypes+1,"constant_pH:protonable"); //ntypes+1 so the atom types start from 1.
   memory->create(typePerProtMol,ntypes+1,"constant_pH:typePerProtMol");
   memory->create(pH1qs,ntypes+1,"constant_pH:pH1qs");
   memory->create(pH2qs,ntypes+1,"constant_pH:pH2qs");



   char line[128];
   if (comm->me == 0)
   {
       if (!pHStructureFile)
           error->all(FLERR,"Error in reading the pH structure file in fix constant_pH");
       fgets(line,sizeof(line),pHStructureFile);
       fgets(line,sizeof(line),pHStructureFile);
       line[strcspn(line,"\n")] = '\0';

       char *token = strtok(line,",");
       pHnTypes = std::stoi(token);
       for (int i = 1; i < ntypes+1; i++)
       {
	   protonable[i] = 0;
	   typePerProtMol[i] = 0;
	   pH1qs[i] = 0.0;
           pH2qs[i] = 0.0;
       }   
       for (int i = 0; i < pHnTypes; i++)
       {
	  if (fgets(line,sizeof(line),pHStructureFile) == nullptr)
	       error->all(FLERR,"Error in reading the pH structure file in fix constant_pH");
	  line[strcspn(line,"\n")] = '\0';
	  token = strtok(line,",");
	  int type = std::stoi(token);
	  protonable[type] = 1;
	  token = strtok(NULL,",");
	  typePerProtMol[type] = std::stoi(token);
	  token = strtok(NULL,",");
	  pH1qs[type] = std::stod(token);
	  token = strtok(NULL,",");
          pH2qs[type] = std::stod(token);
       }
       fclose(pHStructureFile);
       pHStructureFile = nullptr;
   }
   
   MPI_Bcast(protonable,ntypes+1,MPI_INT,0,world);
   MPI_Bcast(typePerProtMol,ntypes+1,MPI_INT,0,world);
   MPI_Bcast(pH1qs,ntypes+1,MPI_DOUBLE,0,world);
   MPI_Bcast(pH2qs,ntypes+1,MPI_DOUBLE,0,world);
}


/* ----------------------------------------------------------------------
   Calculating the number of protonable and HW atoms
   ---------------------------------------------------------------------- */

void FixConstantPH::calculate_num_prot_num_HWs()
{
   double tol = 1e-5;
   int * type = atom->type;
   int nlocal = atom->nlocal;
   int * num_local = new int[2]{0,0};
   int * num_total = new int[2]{0,0};
   
   for (int i = 0; i < nlocal; i++)
   {
      if (type[i] == typeHW)
        num_local[0]++;
      if (protonable[type[i]])
        num_local[1]++;
   }

   MPI_Allreduce(num_local,num_total,2,MPI_INT,MPI_SUM,world);
   num_HWs = num_total[0];
   num_prots = num_total[1];
   
   delete [] num_local;
   delete [] num_total;
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::calculate_dfs()
{
   for (int j = 0; j < n_lambdas; j++) {
	fs[j] = 1.0/(1+exp(-50*(lambdas[j]-0.5)));
        dfs[j] = 50*exp(-50*(lambdas[j]-0.5))*(fs[j]*fs[j]);
   }
}

/* ----------------------------------------------------------------------- */

void FixConstantPH::calculate_dUs()
{
   double U1, U2, U3, U4, U5;
   double dU1, dU2, dU3, dU4, dU5;
   for (int j = 0; j < n_lambdas; j++) {
        U1 = -k*exp(-(lambdas[j]-1-b)*(lambdas[j]-1-b)/(2*a*a));
   	U2 = -k*exp(-(lambdas[j]+b)*(lambdas[j]+b)/(2*a*a));
   	U3 = d*exp(-(lambdas[j]-0.5)*(lambdas[j]-0.5)/(2*s*s));
   	U4 = 0.5*w*(1-erff(r*(lambdas[j]+m)));
   	U5 = 0.5*w*(1+erff(r*(lambdas[j]-1-m)));
   	dU1 = -((lambdas[j]-1-b)/(a*a))*U1;
   	dU2 = -((lambdas[j]+b)/(a*a))*U2;
   	dU3 = -((lambdas[j]-0.5)/(s*s))*U3;
   	dU4 = -0.5*w*r*2*exp(-r*r*(lambdas[j]+m)*(lambdas[j]+m))/sqrt(M_PI);
   	dU5 = 0.5*w*r*2*exp(-r*r*(lambdas[j]-1-m)*(lambdas[j]-1-m))/sqrt(M_PI);

    	Us[j] =  U1 +  U2 +  U3 +  U4 +  U5;
   	dUs[j] = dU1 + dU2 + dU3 + dU4 + dU5;   
   }
}

/* ----------------------------------------------------------------------- */

void FixConstantPH::calculate_dU(const double& _lambda, double& _U, double& _dU)
{
   double U1, U2, U3, U4, U5;
   double dU1, dU2, dU3, dU4, dU5;
   U1 = -k*exp(-(_lambda-1-b)*(_lambda-1-b)/(2*a*a));
   U2 = -k*exp(-(_lambda+b)*(_lambda+b)/(2*a*a));
   U3 = d*exp(-(_lambda-0.5)*(_lambda-0.5)/(2*s*s));
   U4 = 0.5*w*(1-erff(r*(_lambda+m)));
   U5 = 0.5*w*(1+erff(r*(_lambda-1-m)));
   dU1 = -((_lambda-1-b)/(a*a))*U1;
   dU2 = -((_lambda+b)/(a*a))*U2;
   dU3 = -((_lambda-0.5)/(s*s))*U3;
   dU4 = -0.5*w*r*2*exp(-r*r*(_lambda+m)*(_lambda+m))/sqrt(M_PI);
   dU5 = 0.5*w*r*2*exp(-r*r*(_lambda-1-m)*(_lambda-1-m))/sqrt(M_PI);

   _U =  U1 +  U2 +  U3 +  U4 +  U5;
   _dU = dU1 + dU2 + dU3 + dU4 + dU5;   
}


/* ---------------------------------------------------------------------- */

void FixConstantPH::print_Udwp()
{
    double lambda_Udwp, U_Udwp, dU_Udwp;
    int n_points = 100;
	
    lambda_Udwp = -0.5;
    double dlambda_Udwp = 2.0/(double)n_points;

    fprintf(Udwp_fp,"Lambda,U,dU\n");
    for (int i = 0; i <= n_points; i++) {
	calculate_dU(lambda_Udwp,U_Udwp,dU_Udwp);
        fprintf(Udwp_fp,"%f,%f,%f\n",lambda_Udwp,U_Udwp,dU_Udwp);
	lambda_Udwp += dlambda_Udwp;
    }
    fclose(Udwp_fp);
}

/* ----------------------------------------------------------------------
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
  if (q_orig) memory->destroy(q_orig);
  if (f_orig) memory->destroy(f_orig);
  if (peatom_orig) memory->destroy(peatom_orig);
  if (pvatom_orig) memory->destroy(pvatom_orig);
  /* If kspace->force is true these two have been already allocated and 
     here is no need to check it since the lammps destructor first destructs
     the kspace so that "if (force->kspace)" in the destructor for 
     ComputeFEEConstantPH leads to an error!
  */
  
  if (keatom_orig) memory->destroy(keatom_orig);
  if (kvatom_orig) memory->destroy(kvatom_orig);

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
   
void FixConstantPH::modify_qs(double* scales)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;
    double * q = atom->q;


    double * q_changes_local = new double[4]{0.0,0.0,0.0,0.0};
    double * q_changes = new double[4]{0.0,0.0,0.0,0.0};

    // update the charges
    for (int j = 0; j < n_lambdas; j++) {
    	for (int i = 0; i < nlocal; i++)
    	{
	    int molid_i = atom->molecule[i];
            if ((protonable[type[i]] == 1) && (molid_i == molids[j]))
            {
                 double q_init = q_orig[i];
                 q[i] = pH1qs[type[i]] + scales[j] * (pH2qs[type[i]] - pH1qs[type[i]]); // scale == 1 should be for the protonated state
	         q_changes_local[0]++;
	         q_changes_local[1] += (q[i] - q_init);
            }
        }
    }
    
    MPI_Allreduce(q_changes_local,q_changes,2,MPI_DOUBLE,MPI_SUM,world);
    double HW_q_change = -q_changes[1]/static_cast<double>(num_HWs);
    
    for (int i = 0; i < nlocal; i++)
    {
       if (type[i] == typeHW)
       {
	   double q_init = q_orig[i];
	   q[i] = q_init + HW_q_change; //The total charge should be neutral
	   q_changes_local[2]++;
	   q_changes_local[3] += (q[i] - q_init);
       }
    }

    /* The purpose of this part this is just to debug the total charge.
       So, in the final version of the code this part should be 
       commented out!
    */
    /*if (update->ntimestep % nevery == 0) {
    	MPI_Allreduce(q_changes_local,q_changes,4,MPI_DOUBLE,MPI_SUM,world);
    	if (comm->me == 0) error->warning(FLERR,"protonable q change = {}, HW q change = {}, protonable charge change = {}, HW charge change = {}",q_changes[0],q_changes[2],q_changes[1],q_changes[3]);
    }*/
    
    
    //compute_q_total();
    
    delete [] q_changes_local;
    delete [] q_changes;
}


/* ----------------------------------------------------------------------
   modify force and kspace in lammps according
   ---------------------------------------------------------------------- */

void FixConstantPH::update_lmp() {
   int eflag = ENERGY_GLOBAL;
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
   Add forcefield correction term deltaGFF in equation 2 of
   https://pubs.acs.org/doi/full/10.1021/acs.jctc.5b01160
   --------------------------------------------------------------------- */

void FixConstantPH::calculate_GFFs()
{
   for (int j = 0; j < n_lambdas; j++) {
   	int i = 0;
   	while (i < GFF_size && GFF[i][0] < lambdas[j]) i++;

   	if (i == 0)
   	{
      	    error->warning(FLERR,"Warning lambda of {} in Fix constant_pH out of the range, it usually should not happen",lambdas[j]);
            GFF_lambdas[j] = GFF[0][1] + ((GFF[1][1]-GFF[0][1])/(GFF[1][0]-GFF[0][0]))*(lambdas[j] - GFF[0][0]);
        }
        if (i > 0 && i < GFF_size - 1)
            GFF_lambdas[j] = GFF[i-1][1] + ((GFF[i][1]-GFF[i-1][1])/(GFF[i][0]-GFF[i-1][0]))*(lambdas[j] - GFF[i-1][0]);
        if (i == GFF_size - 1)
        {
            error->warning(FLERR,"Warning lambda of {} in Fix constant_pH out of the range, it usually should not happen",lambdas[j]);
            GFF_lambdas[j] = GFF[i][1] + ((GFF[i][1]-GFF[i-1][1])/(GFF[i][0]-GFF[i-1][0]))*(lambdas[j] - GFF[i][0]);
        }
   }
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
   while(fgets(line,sizeof(line), fp) != NULL && i < GFF_size) {
      double _lambda, _GFF;
      line[strcspn(line,"\n")] = 0;
      char * token = strtok(line, ",");
      if (token != NULL) 
	 _lambda = atof(token);
      else 
	 error->one(FLERR,"The GFF correction file in the fix constant_pH is in a wrong format!");
      token = strtok(line, ",");
      if (token != NULL)
	 _GFF = atof(token);
      else
	 error->one(FLERR,"The GFF correction file in the fix constant_pH is in a wrong format!");
      GFF[i][0] = _lambda;
      GFF[i][1] = _GFF;
      i++;   
   }	
   if (fp && (comm->me == 0)) fclose(fp);
}

/* ----------------------------------------------------------------------
   The linear charge interpolation method in Aho et al JCTC 2022
   --------------------------------------------------------------------- */

void FixConstantPH::compute_f_lambda_charge_interpolation()
{
   /* Two different approaches can be used
      either I can go with copying the compute_group_group
      code with factor_lj = 0 or I can use the eng->coul
      I prefer the second one as it is tidier and I guess 
      it should be faster
   */


   int natoms = atom->natoms;
   double * energy_local = new double[n_lambdas];
   double * energy = new double[n_lambdas];
   double * n_lambda_atoms = new double[n_lambdas];

   for (int i = 0; i < n_lambdas; i++) {
      for (int j = 0; j < n_lambda_atoms[i]; j++) {
	  //double delta_q = q_prot[j] - q_deprot[j];
	  // I need to figure out how to identify those atoms
      }
      for (int k = 0; k < n_lambdas; k++) {
	  if (k == i) continue;
	  for (int l = 0; l < n_lambda_atoms[k]; l++) {
	      //double q = (1-lambdas[k])*q_prot[l] + lambdas[k] * q_deprot[l];
              // Double check if the q_prot and q_deprot are in the right place
	      // how should I identify those atoms
	  }
      }
      energy_local[i] = 0.0;
      if (force->pair) energy_local[i] += force->pair->eng_coul;
      // You need to add the kspace contribution too
   }

   MPI_Allreduce(&energy_local, &energy, n_lambdas,MPI_DOUBLE,MPI_SUM,world);
   for (int i = 0; i < n_lambdas; i++) { 
      double force_i = energy[i] / static_cast<double> (natoms); // convert to kcal/mol
      a_lambdas[i] = 4.184*0.0001*force_i / m_lambdas[i]; 
   }
   
   delete [] energy;
   delete [] energy_local;
   delete [] n_lambda_atoms;     
}

/* --------------------------------------------------------------------- */

void FixConstantPH::initialize_v_lambda(const double _T_lambda)
{
    std::mt19937 rng;
    std::normal_distribution<double> distribution(0.0, 1.0);
    double kT = force->boltz * _T_lambda;
    double ke_lambdas = 0.0;
    double ke_lambdas_target = 0.5*n_lambdas*kT; // Not sure about this part.
    for (int j = 0; j < n_lambdas; j++) {
	double stddev = std::sqrt(kT/m_lambdas[j]);
	v_lambdas[j] = distribution(rng);
	ke_lambdas += 0.5*m_lambdas[j]*v_lambdas[j]*v_lambdas[j];
	v_lambdas[j]*= std::sqrt(4184/10.0)/1000.0; // A/fs
    }
    double scaling_factor = std::sqrt(ke_lambdas_target/ke_lambdas);

    for (int j = 0; j < n_lambdas; j++)
	v_lambdas[j] *= scaling_factor;

    double v_cm = 0.0;
    for (int j = 0; j < n_lambdas; j++)
	v_cm += v_lambdas[j];
    v_cm /= static_cast<double>(v_cm);
    for (int j = 0; j < n_lambdas; j++)
	v_lambdas[j] -= v_cm;
}

/* --------------------------------------------------------------------- */

void FixConstantPH::calculate_T_lambda()
{
    T_lambda = 0.0;
    for (int j = 0; j < n_lambdas; j++)
	T_lambda += 0.5*m_lambdas[j]*v_lambdas[j]*v_lambdas[j]*1e7 / (4184*0.0019872041);
}

   
/* --------------------------------------------------------------------- */

void FixConstantPH::compute_q_total()
{
   double * q = atom->q;
   double nlocal = atom->nlocal;
   double q_local = 0.0;
   double tolerance = 0.000001; //0.001;

   for (int i = 0; i <nlocal; i++)
       q_local += q[i];

    MPI_Allreduce(&q_local,&q_total,1,MPI_DOUBLE,MPI_SUM,world);

    if ((q_total >= tolerance || q_total <= -tolerance) && comm->me == 0)
    	error->warning(FLERR,"q_total in fix constant-pH is non-zero: {} from {}",q_total,comm->me);
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
   implementing a method to access the HA, HB, lambda, v_lambda and 
   a_lambda values 
   ---------------------------------------------------------------------- */
   
double FixConstantPH::compute_array(int i, int j)
{ 
   double kj2kcal = 0.239006;
   double kT = force->boltz * T;
   switch(i)
   {
      case 0:
        return HAs[j];
      case 1:
        return HBs[j];
      case 2:
        return dfs[j]*kT*log(10)*(pK-pH);
      case 3:
        return kj2kcal*dUs[j];
      case 4:
        return GFF_lambdas[j];
      case 5:
        return lambdas[j];
      case 6:
        return v_lambdas[j];
      case 7:
        return a_lambdas[j];
      case 8:
        calculate_T_lambda();
        return T_lambda;
      case 9:
        return H_lambdas[j];
   }
   return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array --> Needs to be updated at the end
   ---------------------------------------------------------------------- */

double FixConstantPH::memory_usage()
{
  int nmax = atom->nmax;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  double GFF_bytes = 2.0 * GFF_size * sizeof(double);
  double epsilon_init_bytes = (double) (ntypes + 1) * (double) (ntypes + 1) * sizeof(double);
  double q_orig_bytes = (double) nlocal * sizeof(double);
  double f_orig_bytes = (double) nlocal * 3.0 * sizeof(double);
  double peatom_orig_bytes = (double) nlocal * sizeof(double);
  double pvatom_orig_bytes = (double) nlocal * 6.0 * sizeof(double);
  double keatom_orig_bytes = (double) nlocal * sizeof(double);
  double kvatom_orig_bytes = (double) nlocal * 6.0 * sizeof(double);
  double bytes = GFF_bytes + epsilon_init_bytes + \
	         q_orig_bytes + f_orig_bytes + peatom_orig_bytes + \
                 pvatom_orig_bytes + keatom_orig_bytes + kvatom_orig_bytes;
  return bytes;
}
