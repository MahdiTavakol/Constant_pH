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
/* ---v0.03.02----- */

#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif

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
       lambdas(nullptr), v_lambdas(nullptr), a_lambdas(nullptr), m_lambdas(nullptr), 
       protonable(nullptr), typePerProtMol(nullptr), pH1qs(nullptr), pH2qs(nullptr), q_orig(nullptr), f_orig(nullptr),
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
  
  pstyle = utils::strdup(arg[9]);
  


  qHs = 0.0;
  qHWs = 0.0; //0.278;

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
  
   etha_lambda = 0.0;
   Q_lambda = 100.0;
  
  
   vector_flag = 1;
   size_vector = 10;

}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   #ifdef DEBUG
       std::cout << "Releasing epsilon_init" << std::endl;
   #endif
   memory->destroy(epsilon_init);
   #ifdef DEBUG
       std::cout << "Releasing GFF" << std::endl;
   #endif
   if (GFF != nullptr) memory->destroy(GFF);

   #ifdef DEBUG
       std::cout << "Releasing pparam1" << std::endl;
   #endif
   delete [] pparam1;
   #ifdef DEBUG
       std::cout << "Releasing pstyle" << std::endl;
   #endif
   delete [] pstyle;


   memory->destroy(pH1qs);
   memory->destroy(pH2qs);
   memory->destroy(typePerProtMol);
   memory->destroy(protonable);


   if (lambdas)   delete [] lambdas;
   if (v_lambdas) delete [] v_lambdas;
   if (a_lambdas) delete [] a_lambdas;
   if (m_lambdas) delete [] m_lambdas;

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
   
   pair_params["lj/charmm/coul/charmm"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/gpu"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/intel"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/kk"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/omp"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/implicit"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/implicit/kk"] = "epsilon";
   pair_params["lj/charmm/coul/charmm/implicit/omp"] = "epsilon";
   pair_params["lj/charmm/coul/long"] = "epsilon";
   pair_params["lj/charmm/coul/long/gpu"] = "epsilon";
   pair_params["lj/charmm/coul/long/intel"] = "epsilon";
   pair_params["lj/charmm/coul/long/kk"] = "epsilon";
   pair_params["lj/charmm/coul/long/opt"] = "epsilon";
   pair_params["lj/charmm/coul/long/omp"] = "epsilon";
   pair_params["lj/charmm/coul/msm"] = "epsilon";
   pair_params["lj/charmm/coul/msm/omp"] = "epsilon";
   pair_params["lj/charmmfsw/coul/charmmfsh"] = "epsilon";
   pair_params["lj/charmmfsw/coul/long"] = "epsilon";
   pair_params["lj/charmmfsw/coul/long/kk"] = "epsilon";

   if (pair_params.find(pstyle) == pair_params.end())
      error->all(FLERR,"The pair style {} is not currently supported in fix constant_pH",pstyle);
   
   pparam1 = new char[pair_params[pstyle].length()+1];
   std::strcpy(pparam1,pair_params[pstyle].c_str());
   
   
   
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::setup(int /*vflag*/)
{
   // default values from Donnini, Ullmann, J Chem Theory Comput 2016 - Table S2
    w = 200; //50.0;
    s = 0.3; //0.3;
    h = 4; //0; //10.0;
    k = 2.553; //0; //4.417; //6.267;
    a = 0.034041; //0.04764; //0.04208; //0.05130;
    b = 0.005238; //0.002957; //0.001411;
    r = 16.458; //21.428;
    m = 0.1507; //0.1078;
    d = 2.0; //3.5; //5.0;
    // m_lambda = 20u taken from https://www.mpinat.mpg.de/627830/usage
    m_lambda = 20;
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

    // The limits for these two loops are correct.
    for (int i = 1; i < ntypes+1; i++)
        for (int j = i; j < ntypes+1; j++)
             epsilon_init[i][j] = epsilon[i][j];

	
    GFF_lambda = 0.0;
    GFF = nullptr;
    if (GFF_flag)
	init_GFF();

    if (print_Udwp_flag)
	print_Udwp();


    // Reading the structure of protonable states before and after protonation.
    read_pH_structure_files();


    // Checking if we have enough hydronium ions to neutralize the system
    calculate_num_prot_num_HWs();
	
    fixgpu = modify->get_fix_by_id("package_gpu");


    // This would not work in the initialize section as the m_lambda has not been set yet!
    initialize_v_lambda(310.00);



    n_lambdas = 1; // Just for now
    lambdas   = new double[n_lambdas];
    v_lambdas = new double[n_lambdas];
    a_lambdas = new double[n_lambdas];
    m_lambdas = new double[n_lambdas];

    lambdas[0] = 0.0;
    v_lambdas[0] = 0.0;
    a_lambdas[0] = 0.0;
    m_lambdas[0] = 20.0;

	
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
   calculate_df();
   calculate_dU();
   update_a_lambda();
   compute_Hs<1>();
   compute_q_total();
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
      backup_restore_qfev<1>();      // backup charge, force, energy, virial array values
      modify_q(0.0); //should define a change_parameters(const int);
      update_lmp(); // update the lammps force and virial values
      HA = compute_epair(); 
      backup_restore_qfev<-1>();        // restore charge, force, energy, virial array values
      modify_q(1.0); //should define a change_parameters(const double);
      update_lmp();
      HB = compute_epair();           // HB is for the protonated state with lambda==1 
      backup_restore_qfev<-1>();      // restore charge, force, energy, virial array values
   }
   if (stage == 1)
   {
      modify_q(lambda); //should define a change_parameters(const double);
      //update_lmp(); This update_lmp() might not work here since I am not sure about the correct values for the eflag and vflag variables... Anyway, the epsilon and charge values have been updated according to the pH value and lammps will do the rest
   }
}

/* ----------------------------------------------------------------------
   returns the number of the lambda parameters
  ----------------------------------------------------------------------- */

void FixConstantPH::return_nparams(int& _n_params) const
{
    _n_params = this->n_params;
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

    // just for now
    _x_lambdas[0] = lambda;
    _v_lambdas[0] = v_lambda;
    _a_lambdas[0] = a_lambda;
    _m_lambdas[0] = m_lambda;
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

    // just for now
    lambda = _x_lambdas[0];
    v_lambda = _v_lambdas[0];
    a_lambda = _a_lambdas[0];
    m_lambda = _m_lambdas[0];
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

void FixConstantPH::calculate_df()
{
   f = 1.0/(1+exp(-50*(lambda-0.5)));
   df = 50*exp(-50*(lambda-0.5))*(f*f);
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
   
void FixConstantPH::modify_q(const double& scale)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;
    double * q = atom->q;


    double * q_changes_local = new double[4]{0.0,0.0,0.0,0.0};
    double * q_changes = new double[4]{0.0,0.0,0.0,0.0};

    // update the charges
    for (int i = 0; i < nlocal; i++)
    {
       if (protonable[type[i]] == 1)
       {
           double q_init = q_orig[i];
           q[i] = pH1qs[type[i]] + scale * (pH2qs[type[i]] - pH1qs[type[i]]); // scale == 1 should be for the protonated state
	   q_changes_local[0]++;
	   q_changes_local[1] += (q[i] - q_init);
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

/* ---------------------------------------------------------------------- */

void FixConstantPH::update_a_lambda()
{
   if (GFF_flag) calculate_GFF();
   double NA = 6.022*1e23;
   double kj2kcal = 0.239006;
   double kT = force->boltz * T;

   //df = 1.0;
   //f = 1.0;
   double  f_lambda = -(HB-HA - df*kT*log(10)*(pK-pH) + kj2kcal*dU - GFF_lambda);

   this->a_lambda = 4.184*0.0001*f_lambda / m_lambda;
   /*#ifdef DEBUG
	std::cout << "The a_lambda and f_lambda are :" << a_lambda << "," << f_lambda << std::endl;
   #endif*/

   double  H_lambda = (1-lambda)*HA + lambda*HB - f*kT*log(10*(pK-pH)) + kj2kcal*U + (m_lambda/2.0)*(v_lambda*v_lambda); // This might not be needed. May be I need to tally this into energies.
   // I might need to use the leap-frog integrator and so this function might need to be in other functions than postforce()
	
}

/* --------------------------------------------------------------------- */

void FixConstantPH::initialize_v_lambda(const double _T_lambda)
{
    v_lambda = sqrt(0.0019872041*4184.0 * _T_lambda / (10.0 * m_lambda))/1000.0;
}

/* --------------------------------------------------------------------- */

void FixConstantPH::calculate_T_lambda()
{
    T_lambda = m_lambda * v_lambda * v_lambda*1e7 / (4184*0.0019872041);
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
   
double FixConstantPH::compute_vector(int i)
{ 
   double kj2kcal = 0.239006;
   double kT = force->boltz * T;
   switch(i)
   {
      case 0:
        return HA;
      case 1:
        return HB;
      case 2:
        return df*kT*log(10)*(pK-pH);
      case 3:
        return kj2kcal*dU;
      case 4:
        return GFF_lambda;
      case 5:
        return lambda;
      case 6:
        return v_lambda;
      case 7:
        return a_lambda;
      case 8:
        calculate_T_lambda();
        return T_lambda;
      case 9:
        return a_etha_v_ratio_lambda;
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
