
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
/* ---v0.07.04----- */

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
#include "random_park.h"

#include <cstring>
#include <string>
#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { 
       NONE=0,
       BUFFER=1<<0, 
       ADAPTIVE=1<<1,
       ZEROCHARGE=1<<2,
       CONSTRAIN=1<<3
     };

static constexpr double tol = 1e-5;
/* ---------------------------------------------------------------------- */

FixConstantPH::FixConstantPH(LAMMPS *lmp, int narg, char **arg): Fix(lmp, narg, arg),
       pHStructureFile(nullptr), 
       pH1qs(nullptr), pH2qs(nullptr), typePerProtMol(nullptr), protonable(nullptr),
       HAs(nullptr), HBs(nullptr), Us(nullptr), dUs(nullptr),
       lambdas(nullptr), v_lambdas(nullptr), a_lambdas(nullptr), m_lambdas(nullptr), H_lambdas(nullptr),
       molids(nullptr),
       fs(nullptr), dfs(nullptr),
       fp(nullptr), GFF(nullptr), GFF_lambdas(nullptr),
       Udwp_fp(nullptr), fix_adaptive_protonation_id(nullptr),
       fixgpu(nullptr), q_orig(nullptr), f_orig(nullptr),
       peatom_orig(nullptr), pvatom_orig(nullptr),keatom_orig(nullptr), kvatom_orig(nullptr)
{
  if (narg < 8) utils::missing_cmd_args(FLERR,"fix constant_pH", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Illegal fix constant_pH every value {}", nevery);
  // Reading the file that contains the charges before and after protonation/deprotonation
  if (comm->me == 0) {
      pHStructureFile = fopen(arg[4],"r"); // The command reads the file the type and charge of each atom before and after protonation
      if (pHStructureFile == nullptr)
         error->all(FLERR,"Unable to open the file");
  }
  
  

	
  pK = utils::numeric(FLERR, arg[5], false, lmp);
  pH = utils::numeric(FLERR, arg[6], false, lmp);
  T = utils::numeric(FLERR, arg[7], false, lmp);
  


  qOWs = -0.834;
  qHWs = 0.278;
  
  // The default value for the mu
  mu = 0.0;

  /* Unset all the flags
     it is an important step since
     in C++ it is not guaranteed that
     the default value of an int to be
     zero 
     */
  flags = 0;

  GFF_flag = false;
  print_Udwp_flag = false;
  n_lambdas = 1;
  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "GFF") == 0)
    {
       GFF_flag = true;
       fp = fopen(arg[iarg+1],"r");
       if (fp == nullptr)
         error->one(FLERR, "Cannot find fix constant_pH the GFF correction file {}",arg[iarg+1]);
       iarg = iarg + 2;
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
	iarg+=2;
	memory->create(molids,n_lambdas,"constant_pH:lambdas");
	for (int i = 0; i < n_lambdas; i++) {
	    molids[i] = utils::numeric(FLERR,arg[iarg],false,lmp);
	    iarg++;
	}
    }
    else if (strcmp(arg[iarg],"mu") == 0) {
        if (narg < iarg+2) utils::missing_cmd_args(FLERR,"fix constant_pH", error);
        mu = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        iarg += 2;
        }
    else if (strcmp(arg[iarg],"buffer") == 0) {
	flags |= BUFFER;
	if (narg < iarg+6) utils::missing_cmd_args(FLERR,"fix constant_pH", error);
	N_buff  = utils::numeric(FLERR,arg[iarg+1],false,lmp);
	typeOWs = utils::numeric(FLERR,arg[iarg+2],false,lmp);
	typeHWs = utils::numeric(FLERR,arg[iarg+3],false,lmp);
	if (typeOWs > atom->ntypes) error->all(FLERR,"Illegal fix constant_pH atom type {}",typeOWs);
	if (typeHWs > atom->ntypes) error->all(FLERR,"Illegal fix constant_pH atom type {}",typeHWs);
	qOWs = utils::numeric(FLERR,arg[iarg+4],false,lmp);
	qHWs = utils::numeric(FLERR,arg[iarg+5],false,lmp);
	iarg += 6;
    }
    else if (strcmp(arg[iarg],"Fix_adaptive_protonation") == 0) {
	flags |= ADAPTIVE;
	fix_adaptive_protonation_id = utils::strdup(arg[iarg+1]);
	nevery_fix_adaptive = utils::numeric(FLERR,arg[iarg+2],false,lmp);
	iarg+=3;
    }
    else if (strcmp(arg[iarg],"constrain") == 0) {
        flags |= CONSTRAIN;
        iarg++;
    }
    else if (strcmp(arg[iarg],"zero_total_charge") == 0) {
        flags |= ZEROCHARGE;
        iarg++;
    }
    else
       error->all(FLERR, "Unknown fix constant_pH keyword: {}", arg[iarg]);
   }
  
   fixgpu = nullptr;
  
  
  
   array_flag = 1;
   size_array_rows = 12;
   size_array_cols = n_lambdas+1;

}

/* ---------------------------------------------------------------------- */

FixConstantPH::~FixConstantPH()
{
   // Closing files if they are open. 
   if (pHStructureFile && comm->me == 0) fclose(pHStructureFile); // I have already closed it.. This is here just to assure it is closed
   if (fp && (comm->me == 0)) fclose(fp);
   if (Udwp_fp && (comm->me == 0)) fclose(Udwp_fp); // We should never reach that point as this file is writting just at the setup stage and then it will be closed

   // deallocate char* variables
   if (fix_adaptive_protonation_id) delete [] fix_adaptive_protonation_id; // Since it is allocated with lmp->utils->strdup, it must be deallocated with delete [] 

   // deallocating the variables whose size depend on the ntypes and as the ntypes does not change during the simulation there is no need for reallocation of them
   if (pH1qs) memory->destroy(pH1qs);
   if (pH2qs) memory->destroy(pH2qs);
   if (typePerProtMol) memory->destroy(typePerProtMol);
   if (protonable) memory->destroy(protonable);
   if (GFF) memory->destroy(GFF);

   // deallocate the memories with size dependent on the n_lambda	
   delete_lambdas();

   // deallocate memories whose size is dependent on natoms
   deallocate_storage();   
}

/* ----------------------------------------------------------------------

   ---------------------------------------------------------------------- */

int FixConstantPH::setmask()
{
   int mask = 0;
   mask |= INITIAL_INTEGRATE; // Calculates the a_lambda
   mask |= POST_FORCE; // Updates the a_lambda before the second step of the velocity verlet
   return mask;	
}

/* ----------------------------------------------------------------------
   Setup
   ---------------------------------------------------------------------- */
   
void FixConstantPH::init()
{
   if (flags & ADAPTIVE) {
      fix_adaptive_protonation = static_cast<FixAdaptiveProtonation*>(modify->get_fix_by_id(fix_adaptive_protonation_id));
      fix_adaptive_protonation->get_n_protonable(n_lambdas);
   }
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::setup(int /*vflag*/)
{
   // default values from Donnini, Ullmann, J Chem Theory Comput 2016 - Table S2
    w = 200;
    s = 0.3;
    h = 7.0;
    k = 2.553; //4.417; //2.553;
    a = 0.03401; //0.04208; //0.03401;
    b = 0.005238; //0.002957; //0.005238;
    r = 16.458; 
    m = 0.1507;
    d = 2.0; //5.00; //2.0; //The height of the barrier is 2*d

    // default values for the buffer potential with h = 0 from Donnin J Chem Theory Comput 2016 - Table S2
    w_buff = 200;
    s_buff = 0.3;
    h_buff = 0.0;
    k_buff = 0.0;
    a_buff = 0.04764;
    b_buff = -0.09706;
    r_buff = 16.458;
    m_buff = 0.1507;
    d_buff = 0.0;
	
    // Reading the structure of protonable states before and after protonation.
    read_pH_structure_files();


    // Checking if we have correct number of hydronium ions
    if (flags & BUFFER) check_num_OWs_HWs();
	
    fixgpu = modify->get_fix_by_id("package_gpu");
    
    
    /* As it is hypothesized that the initial values for 
     * lambdas are zero the initial value for the lambda_buff
     * should be one so there is enough protons to be exchanged
     * between the lambdas and buffer due to the contraint on
     * the lambas[0] + ... + lambdas[n] + lambda_buff
     */
    
    if (flags & BUFFER) {
        lambda_buff = 1.0;
        v_lambda_buff = 0.0;
        m_lambda_buff = 20.0;
        
        modify_q_buff(lambda_buff);
        compute_q_total();
    }


    set_lambdas();
       
    if (GFF_flag)
	init_GFF();

    if (print_Udwp_flag)
	print_Udwp();


    nmax = atom->nmax;
    allocate_storage();
    
}

/* ----------------------------------------------------------------------
   This part calculates the acceleration of the lambdas parameter
   which is obtained from the force acting on it
   ---------------------------------------------------------------------- */

void FixConstantPH::initial_integrate(int /*vflag*/)
{
   if (flags & ADAPTIVE) {
      if (!(update->ntimestep % nevery_fix_adaptive))
      {
         int n_protonable;
         fix_adaptive_protonation->get_n_protonable(n_protonable);
         if (n_protonable != this->n_lambdas) {
            this->n_lambdas = n_protonable;
            delete_lambdas();
            set_lambdas();
            fix_adaptive_protonation->get_protonable_molids(molids);
         }       
      }
   }
   	
   	
   calculate_dfs();
   calculate_dUs();
   update_a_lambda();
}

/* ----------------------------------------------------------------------
   The second step of the integration 
   ----------------------------------------------------------------------  */
   
void FixConstantPH::post_force(int /*vflag*/)
{
   calculate_dfs();
   calculate_dUs();
   update_a_lambda();
}

/* ----------------------------------------------------------------------
   This function deallocates the storage for memories whose sizes are 
   dependent on the n_lambdas
   ----------------------------------------------------------------------  */

void FixConstantPH::delete_lambdas()
{
   if (lambdas) memory->destroy(lambdas);
   if (v_lambdas) memory->destroy(v_lambdas);
   if (a_lambdas) memory->destroy(a_lambdas);
   if (m_lambdas) memory->destroy(m_lambdas);
   if (H_lambdas) memory->destroy(H_lambdas);

   if (HAs) memory->destroy(HAs);
   if (HBs) memory->destroy(HBs);
   if (GFF_lambdas) memory->destroy(GFF_lambdas);
   if (fs) memory->destroy(fs);
   if (dfs) memory->destroy(dfs);
   if (Us) memory->destroy(Us);
   if (dUs) memory->destroy(dUs);
   
   if (lambdas_j) memory->destroy(lambdas_j);

   if (molids) memory->destroy(molids);
}

/* ----------------------------------------------------------------------
   This function allocates the storage for memories whose sizes are 
   dependent on the n_lambdas
   ----------------------------------------------------------------------  */

void FixConstantPH::set_lambdas() {
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
   
   memory->create(lambdas_j,n_lambdas,"constant_pH:lambdas_j");

   if (flags & ADAPTIVE) memory->create(molids,n_lambdas,"constant_pH:molids");

   for (int i = 0; i < n_lambdas; i++) {
      lambdas[i] = 0.0;
      v_lambdas[i] = 0.0;
      a_lambdas[i] = 0.0;
      m_lambdas[i] = 20.0; // m_lambda == 20.0u taken from https://www.mpinat.mpg.de/627830/usage
      H_lambdas[i] = 0.0;
      GFF_lambdas[i] = 0.0;
      if (flags & ADAPTIVE)
         molids[i] = 0.0;
   } 

   // This would not work in the initialize section as the m_lambda has not been set yet!
   initialize_v_lambda(this->T);
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::update_a_lambda()
{
   if (GFF_flag) calculate_GFFs();
   double mvv2e = force->mvv2e;
   double kj2kcal = 0.239006;
   double kT = force->boltz * T;

   //df = 1.0;
   //f = 1.0;

   for (int i = 0; i < n_lambdas; i++) {
	double  f_lambda = -(-dfs[i]*kT*log(10)*(pK-pH) + kj2kcal*dUs[i] - GFF_lambdas[i]); // The df sign should be positive if the lambda = 0 is for the protonated state 
	this->a_lambdas[i] = f_lambda /m_lambdas[i]; // 4.184*0.0001*f_lambda / m_lambda;
	// I am not sure about the sign of the f*kT*log(10)*(pK-pH)
        this->H_lambdas[i] = - fs[i]*kT*log(10)*(pK-pH) + kj2kcal*Us[i] + (m_lambdas[i]/2.0)*(v_lambdas[i]*v_lambdas[i])*mvv2e; // This might not be needed. May be I need to tally this into energies.
        // I might need to use the leap-frog integrator and so this function might need to be in other functions than postforce()
   }

   if (flags & BUFFER) {
	double f_lambda_buff = -(kj2kcal*dU_buff);
	this->a_lambda_buff = f_lambda_buff / m_lambda_buff; // the fix_nh_constant_pH itself takes care of units
	this->H_lambda_buff = kj2kcal*U_buff + N_buff*(m_lambda_buff/2.0)*(v_lambda_buff*v_lambda_buff)*mvv2e;
   }
}

/* ----------------------------------------------------------------------- 
    This function is called by compute_GFF for thermodynamic integration 
    of the GFF value.
   ----------------------------------------------------------------------- */

void FixConstantPH::calculate_H_once()
{
   calculate_dfs();
   calculate_dUs();
   update_a_lambda();
}
	
/* ----------------------------------------------------------------------- */

void FixConstantPH::compute_Hs()
{
   if (nmax < atom->nmax)
   {
      nmax = atom->nmax;
      allocate_storage();
      deallocate_storage();
   }
   
   backup_restore_qfev<1>();
   // computing the HA and HB for each lambda
   for (int j = 0; j < n_lambdas; j++) {
      std::fill(lambdas_j,lambdas_j+n_lambdas,0.0);
      double lambda_j = 0.0;
      modify_qs(lambda_j,j);
      update_lmp();
      HAs[j] = compute_epair();
      backup_restore_qfev<-1>();
      lambda_j = 1.0;
      modify_qs(lambda_j,j);
      update_lmp();
      HBs[j] = compute_epair();
      backup_restore_qfev<-1>();
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
    This function sets the value of qs based on the value of lambdas()
    and lambda_buffs() calculated by the nve_x() fuction of the 
    fix_nh_constant_pH. Accordingly, this function should be called from
    the last line of the nve_x() fuction.
   --------------------------------------------------------------------- */

void FixConstantPH::reset_qs()
{
    modify_qs(lambdas);
    
    if (flags & BUFFER) modify_q_buff(lambda_buff);
        
    /* This should be here just for debugging
       since it used MPI_Allreduce to calculate
       the total charge it has some overhead not 
       advised in the production run
    */
    /*
       There is no need for this anymore 
       since the q_total = sigma_lambdas * mol_charge_change + N_buff* lambda_buff*buff_charge_change
    */
    if (0) {
        compute_q_total();
    }
}

/* ---------------------------------------------------------------------
    This function returns the H_lambdas
   --------------------------------------------------------------------- */

void FixConstantPH::return_H_lambdas(double* _H_lambdas) const
{
    for (int i = 0; i < n_lambdas; i++)
        _H_lambdas[i] = H_lambdas[i]; 
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

/* ----------------------------------------------------------------------
   returns the buffer parameters
   ---------------------------------------------------------------------- */

void FixConstantPH::return_buff_params(double& _x_lambda_buff, double& _v_lambda_buff, 
                                         double& _a_lambda_buff, double& _m_lambda_buff, int& _N_buff) const
{
    if (!(flags & BUFFER)) {
	error->warning(FLERR,"There is no buffer in the fix constant pH so you should not have reached here");
	return;
    }
    _x_lambda_buff = this->lambda_buff;
    _v_lambda_buff = this->v_lambda_buff;
    _a_lambda_buff = this->a_lambda_buff;
    _m_lambda_buff = this->m_lambda_buff;
    _N_buff = this->N_buff;
}

/* ----------------------------------------------------------------------
   sets the buffer parameters which should be used by the fix_nh_constant_pH
   the N_buff is not here since I do not want other commands to mess up with
   the number of buffers
   ---------------------------------------------------------------------- */

void FixConstantPH::reset_buff_params(const double _x_lambda_buff, const double _v_lambda_buff, 
                                        const double _a_lambda_buff, const double _m_lambda_buff) 
{
    if (!(flags & BUFFER)) {
	error->warning(FLERR,"There is no buffer in the fix constant pH so you should not have reached here");
	return;
    }
    this->lambda_buff = _x_lambda_buff;
    this->v_lambda_buff = _v_lambda_buff;
    this->a_lambda_buff = _a_lambda_buff;
    this->m_lambda_buff = _m_lambda_buff;
}

/* ---------------------------------------------------------------------- */

void FixConstantPH::read_pH_structure_files()
{
   /* File format
    * Comment 
    * pHnTypes
    * type1,  number of type1 atoms in the protonable molecule, qBeforeProtonation, qAfterProtonation
    * ... 
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
   }
   
   pHStructureFile = nullptr;
   
   MPI_Bcast(protonable,ntypes+1,MPI_INT,0,world);
   MPI_Bcast(typePerProtMol,ntypes+1,MPI_INT,0,world);
   MPI_Bcast(pH1qs,ntypes+1,MPI_DOUBLE,0,world);
   MPI_Bcast(pH2qs,ntypes+1,MPI_DOUBLE,0,world);
}


/* ----------------------------------------------------------------------
   Checking the number of Oxygen and hydrogen atoms of hydronium ions in 
   the simulation box.
   ---------------------------------------------------------------------- */

void FixConstantPH::check_num_OWs_HWs()
{
   int * type = atom->type;
   int nlocal = atom->nlocal;
   int * num_local = new int[2]{0,0};
   int * num_total = new int[2]{0,0};
   
   for (int i = 0; i < nlocal; i++)
   {
      if (type[i] == typeHWs)
        num_local[0]++;
      if (type[i] == typeOWs)
        num_local[1]++;
   }

   MPI_Allreduce(num_local,num_total,2,MPI_INT,MPI_SUM,world);
   num_HWs = num_total[0];
   num_OWs = num_total[1];

  if (num_HWs != 3*num_OWs) 
      error->one(FLERR,"Number of HWs in the fix constant pH {} is not three times the number of OWs {}",num_HWs,num_OWs);
  if (num_OWs != N_buff)
      error->one(FLERR,"Wrong number of N_buff in the fix constant pH: {}",N_buff);
   
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
        U1 = -k*exp(-(lambdas[j]-1.0-mu-b)*(lambdas[j]-1.0-mu-b)/(2.0*a*a));
   	U2 = -k*exp(-(lambdas[j]+mu+b)*(lambdas[j]+mu+b)/(2.0*a*a));
   	U3 = d*exp(-(lambdas[j]-0.5)*(lambdas[j]-0.5)/(2.0*s*s));
   	U4 = 0.5*w*(1.0-erff(r*(lambdas[j]+m)));
   	U5 = 0.5*w*(1.0+erff(r*(lambdas[j]-1.0-m)));
   	dU1 = -((lambdas[j]-1.0-b)/(a*a))*U1;
   	dU2 = -((lambdas[j]+b)/(a*a))*U2;
   	dU3 = -((lambdas[j]-0.5)/(s*s))*U3;
   	dU4 = -0.5*w*r*2*exp(-r*r*(lambdas[j]+m)*(lambdas[j]+m))/sqrt(M_PI);
   	dU5 = 0.5*w*r*2*exp(-r*r*(lambdas[j]-1-m)*(lambdas[j]-1.0-m))/sqrt(M_PI);

    	Us[j] =  U1 +  U2 +  U3 +  U4 +  U5;
   	dUs[j] = dU1 + dU2 + dU3 + dU4 + dU5;   
   }

   if (flags & BUFFER) {
	U1 = -k_buff*exp(-(lambda_buff-1.0-b_buff)*(lambda_buff-1.0-b)/(2.0*a_buff*a_buff));
   	U2 = -k_buff*exp(-(lambda_buff+b_buff)*(lambda_buff+b_buff)/(2.0*a_buff*a_buff));
   	U3 = d_buff*exp(-(lambda_buff-0.5)*(lambda_buff-0.5)/(2*s_buff*s_buff));
   	U4 = 0.5*w_buff*(1.0-erff(r_buff*(lambda_buff+m_buff)));
   	U5 = 0.5*w_buff*(1.0+erff(r_buff*(lambda_buff-1.0-m_buff)));
   	dU1 = -((lambda_buff-1.0-b_buff)/(a_buff*a_buff))*U1;
   	dU2 = -((lambda_buff+b_buff)/(a_buff*a_buff))*U2;
   	dU3 = -((lambda_buff-0.5)/(s_buff*s_buff))*U3;
   	dU4 = -0.5*w_buff*r_buff*2*exp(-r_buff*r_buff*(lambda_buff+m_buff)*(lambda_buff+m_buff))/sqrt(M_PI);
   	dU5 = 0.5*w_buff*r_buff*2*exp(-r_buff*r_buff*(lambda_buff-1.0-m_buff)*(lambda_buff-1.0-m_buff))/sqrt(M_PI);

    	U_buff =  U1 +  U2 +  U3 +  U4 +  U5;
   	dU_buff = dU1 + dU2 + dU3 + dU4 + dU5;  	   
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

    if (comm->me == 0) {
        fprintf(Udwp_fp,"Lambda,U,dU\n");
        for (int i = 0; i <= n_points; i++) {
	    calculate_dU(lambda_Udwp,U_Udwp,dU_Udwp);
            fprintf(Udwp_fp,"%f,%f,%f\n",lambda_Udwp,U_Udwp,dU_Udwp);
	    lambda_Udwp += dlambda_Udwp;
        }
        fclose(Udwp_fp);
    }
    Udwp_fp = nullptr;
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
  memory->create(q_orig, nlocal, "constant_pH:q_orig");
  memory->create(f_orig, nlocal, 3, "constant_pH:f_orig");
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

  q_orig = nullptr;
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
   modify just q of one lambda
   -------------------------------------------------------------- */
   
void FixConstantPH::modify_qs(double scale, int j)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;
    double * q = atom->q;


    double * q_changes_local = new double[4]{0.0,0.0,0.0,0.0};
    double * q_changes = new double[4]{0.0,0.0,0.0,0.0};
    
    
    for (int i = 0; i < nlocal; i++) {
        int molid_i = atom->molecule[i];
        if ((protonable[type[i]] == 1) && (molid_i == molids[j]))
        {
            double q_init = q_orig[i];
            q[i] = pH1qs[type[i]] + scale * (pH2qs[type[i]] - pH1qs[type[i]]); // scale == 1 should be for the protonated state
	    q_changes_local[0]++;
	    q_changes_local[1] += (q[i] - q_init);
        }
    }


    /* If the buffer is set the modify_q_buffer modifies the charge of the buffer 
       and the constraint in the fix_nh_constant_pH would constrain the total charge.
       So, nothing lefts to do here! */
    if (!(flags & BUFFER) || (flags & ZEROCHARGE)) {
    	MPI_Allreduce(q_changes_local,q_changes,2,MPI_DOUBLE,MPI_SUM,world);
	double HW_q_change = -q_changes[1]/static_cast<double>(num_HWs);
	
	

	for (int i = 0; i < nlocal; i++) {
            if (type[i] == typeHWs) {
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
        }
        compute_q_total();*/
    }
    
    delete [] q_changes_local;
    delete [] q_changes;
}

/* --------------------------------------------------------------
   modify the q of the lambdas
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



    /* If the buffer is set the modify_q_buffer modifies the charge of the buffer 
       and the constraint in the fix_nh_constant_pH would constrain the total charge.
       So, nothing lefts to do here! */
    if (!(flags & BUFFER) || (flags & ZEROCHARGE)) {
    	MPI_Allreduce(q_changes_local,q_changes,2,MPI_DOUBLE,MPI_SUM,world);
	double HW_q_change = -q_changes[1]/static_cast<double>(num_HWs);
	
	

	for (int i = 0; i < nlocal; i++) {
            if (type[i] == typeHWs) {
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
        }
        compute_q_total();*/
    }
    
    delete [] q_changes_local;
    delete [] q_changes;
}

/* --------------------------------------------------------------
   modify the q of the buffer
   -------------------------------------------------------------- */
   
void FixConstantPH::modify_q_buff(const double _scale)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;
    double * q = atom->q;


    // update the charges
    
    for (int i = 0; i < nlocal; i++)
    {
        if (type[i] == typeHWs) {
	    q[i] = (_scale-qOWs) / 3.0;
        } else if (type[i] == typeOWs) {
            q[i] = qOWs; // Just to assure if the charge of Oxygen atoms of the hydronium ions are correct!
        }
    }
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
   fp = nullptr;
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

   delete [] energy_local;
   delete [] energy;
   delete [] n_lambda_atoms;     
}

/* --------------------------------------------------------------------- */

void FixConstantPH::initialize_v_lambda(const double _T_lambda)
{
    RanPark *random = nullptr;
    double seed = 1234579;
    random = new RanPark(lmp,seed);

    
    for (int j = 0; j < n_lambdas; j++) {
	v_lambdas[j] = random->gaussian()/std::sqrt(m_lambdas[j]);
    }
    if (flags & BUFFER) {
        v_lambda_buff = random->gaussian()/std::sqrt(m_lambda_buff);
    }
    
    
    this->calculate_T_lambda();

	
    double scaling_factor = std::sqrt(_T_lambda/T_lambda);

    for (int j = 0; j < n_lambdas; j++)
	v_lambdas[j] *= scaling_factor;
    if (flags & BUFFER)
        v_lambda_buff *= scaling_factor;
        
    
    this->calculate_T_lambda();    
    scaling_factor = std::sqrt(_T_lambda/T_lambda);
    
    

    double v_cm = 0.0;
    for (int j = 0; j < n_lambdas; j++)
	v_cm += v_lambdas[j];
    if (flags & BUFFER)
        v_cm += N_buff*v_lambda_buff;
    double n_cm = static_cast<double>(n_lambdas);
    if (flags & BUFFER) n_cm += static_cast<double>(N_buff);
    v_cm /= static_cast<double>(n_cm);
    for (int j = 0; j < n_lambdas; j++)
	v_lambdas[j] -= v_cm;
    if (flags & BUFFER) 
        v_lambda_buff -= v_cm;
	
	
	
    delete random;
}

/* --------------------------------------------------------------------- */

void FixConstantPH::calculate_T_lambda()
{
    double KE_lambda = 0.0;
    double KE_lambda_buff = 0.0;
    double k = force->boltz;
    double mvv2e = force->mvv2e;
    
    double Nf = static_cast<double>(n_lambdas);
    if (flags & BUFFER)
    	Nf += static_cast<double>(N_buff);
    if (flags & CONSTRAIN)
    	Nf -= 1.0;
    	
    for (int j = 0; j < n_lambdas; j++)
        KE_lambda += 0.5*m_lambdas[j]*v_lambdas[j]*v_lambdas[j]*mvv2e;
    if (flags & BUFFER)
        KE_lambda_buff += 0.5*N_buff*m_lambda_buff*v_lambda_buff*v_lambda_buff*mvv2e;
    
    T_lambda = 2.0*(KE_lambda+KE_lambda_buff) / (Nf * k);
    
}

   
/* --------------------------------------------------------------------- */

double FixConstantPH::compute_q_total()
{
   double * q = atom->q;
   double nlocal = atom->nlocal;
   double q_local = 0.0;
   double tolerance = 1e-6; //0.001;

   for (int i = 0; i <nlocal; i++)
      q_local += q[i];

   MPI_Allreduce(&q_local,&q_total,1,MPI_DOUBLE,MPI_SUM,world);
      
   
   if (std::abs(q_total) > tolerance && comm->me == 0)
      error->warning(FLERR,"q_total in fix constant-pH is non-zero: {} from {}",q_total,comm->me);
      
   return q_total;
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
        // 1
        if (j < n_lambdas)
           return HAs[j];
        else if (j == n_lambdas)
           return N_buff*HA_buff;
        else
           return -1.0;
      case 1:
        // 2
        if (j < n_lambdas)
           return HBs[j];
        else if (j == n_lambdas)
           return N_buff*HB_buff;
        else
           return -1.0;
      case 2:
        // 3
        if (j < n_lambdas)
           return dfs[j]*kT*log(10)*(pK-pH);
        else if (j == n_lambdas)
           return 0.0;
        else
           return -1.0;
      case 3:
        // 4
        if (j < n_lambdas)
           return kj2kcal*dUs[j];
        else if (j == n_lambdas)
           return dU_buff;
        else
           return -1.0;
      case 4:
        // 5
        if (j < n_lambdas)
           return GFF_lambdas[j];
        else if (j == n_lambdas)
           return 0.0;
        else
           return -1.0;
      case 5:
        // 6
        if (j < n_lambdas)
           return lambdas[j];
        else if (j == n_lambdas)
           return lambda_buff;
        else
           return -1.0;
      case 6:
        // 7
        if (j < n_lambdas)
           return v_lambdas[j];
        else if (j == n_lambdas)
           return v_lambda_buff;
        else
           return -1.0;
      case 7:
        // 8
        if (j < n_lambdas)
           return a_lambdas[j];
        else if (j == n_lambdas)
           return a_lambda_buff;
        else
           return -1.0;
      case 8:
        // 9
        calculate_T_lambda();
        return T_lambda;
      case 9:
        // 10
        if (j < n_lambdas)
           return H_lambdas[j];
        else if (j == n_lambdas)
           return H_lambda_buff;
        else
           return -1.0;
      case 10:
        // 11
        compute_q_total();
        return q_total;
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
