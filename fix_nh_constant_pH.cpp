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

/* ----------------------------------------------------------------------
   Contributing authors: Mark Stevens (SNL), Aidan Thompson (SNL)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Constant pH support added by: Mahdi Tavakol (Oxford)
   v0.05.19
------------------------------------------------------------------------- */
#include "fix_constant_pH.h"
#include "fix_nh_constant_pH.h"


#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "force.h"
#include "group.h"
#include "irregular.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <random>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double DELTAFLIP = 0.1;
static constexpr double TILTMAX = 1.5;
static constexpr double EPSILON = 1.0e-6;

enum{NOBIAS,BIAS};
enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

// enums for the lambda integration
enum {LAMBDA_ANDERSEN,LAMBDA_BOSSI,LAMBDA_NOSEHOOVER};
enum { 
       NONE_LAMBDA=0,
       BUFFER=1<<0, 
       CONSTRAIN=1<<1
     };

/* ----------------------------------------------------------------------
   NVT,NPH,NPT integrators for improved Nose-Hoover equations of motion
 ---------------------------------------------------------------------- */

FixNHConstantPH::FixNHConstantPH(LAMMPS *lmp, int narg, char **arg) :
    FixNH(lmp, narg, arg), 
    fix_constant_pH(nullptr), fix_constant_pH_id(nullptr), 
    x_lambdas(nullptr), v_lambdas(nullptr), a_lambdas(nullptr), m_lambdas(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);
  
  lambda_thermostat_type = NONE_LAMBDA;
  lambda_integration_flags = 0;
  

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"fix_constant_pH_id") == 0) {
       lambda_thermostat_type = NONE_LAMBDA;
       lambda_integration_flags = 0;
       fix_constant_pH_id = utils::strdup(arg[iarg+1]);
       iarg += 2;
    } else if (strcmp(arg[iarg],"lambda_andersen") == 0) {
       lambda_thermostat_type = LAMBDA_ANDERSEN;
       t_andersen = utils::numeric(FLERR,arg[iarg+1],false,lmp);
       iarg+=2;
    } else if (strcmp(arg[iarg],"lambda_bussi") == 0) {
       lambda_thermostat_type = LAMBDA_BOSSI;
       iarg++;
    } else if (strcmp(arg[iarg],"lambda_nose-hoover") == 0) {
       lambda_thermostat_type = LAMBDA_NOSEHOOVER;
       Q_lambda_nose_hoover = utils::numeric(FLERR,arg[iarg+1],false,lmp);
       iarg+=2;
    } else if (strcmp(arg[iarg],"buffer") == 0) {
       lambda_integration_flags |= BUFFER;
       iarg++;
    } else if (strcmp(arg[iarg],"constrain_total_charge") == 0) {
       if (!(lambda_integration_flags & BUFFER))
          error->one(FLERR,"Constrain total charge in absence of a buffer is not supported yet!");
       lambda_integration_flags |= CONSTRAIN;
       mols_charge_change = utils::numeric(FLERR,arg[iarg+1],false,lmp);
       buff_charge_change = utils::numeric(FLERR,arg[iarg+2],false,lmp);
       total_charge = utils::numeric(FLERR,arg[iarg+3],false,lmp);
       iarg += 4;
    } else {
       // skip to next argument; argument check for unknown keywords is done in FixNH
       ++iarg;
    }
  }

  if (fix_constant_pH_id == nullptr) error->all(FLERR, "Invalid fix_nh constant_pH");

}

/* ---------------------------------------------------------------------- */

FixNHConstantPH::~FixNHConstantPH()
{
  FixNH::~FixNH();
  if (fix_constant_pH_id) delete [] fix_constant_pH_id;
  if (x_lambdas) memory->destroy(x_lambdas);
  if (v_lambdas) memory->destroy(v_lambdas);
  if (a_lambdas) memory->destroy(a_lambdas);
  if (m_lambdas) memory->destroy(m_lambdas); 
}

/* ---------------------------------------------------------------------- */

void FixNHConstantPH::init()
{
  FixNH::init();

  
  fix_constant_pH = static_cast<FixConstantPH*>(modify->get_fix_by_id(fix_constant_pH_id));
  fix_constant_pH->return_nparams(n_lambdas);

   
  memory->create(x_lambdas,n_lambdas,"nh_constant_pH:x_lambdas");
  memory->create(v_lambdas,n_lambdas,"nh_constant_pH:v_lambdas");
  memory->create(a_lambdas,n_lambdas,"nh_constant_pH:a_lambdas");
  memory->create(m_lambdas,n_lambdas,"nh_constant_pH:m_lambdas");

  std::srand(static_cast<unsigned int>(std::time(nullptr)));

  zeta_bussi = 0.0;
  zeta_nose_hoover = 0.0;
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities
-----------------------------------------------------------------------*/

void FixNHConstantPH::nve_v()
{
  FixNH::nve_v();
  
  bigint ntimestep = update->ntimestep;

  fix_constant_pH->return_nparams(n_lambdas);
  fix_constant_pH->return_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);

  for (int i = 0; i < n_lambdas; i++)
     v_lambdas[i] += dtf * a_lambdas[i];
  fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);

  if (lambda_integration_flags & BUFFER) {
     fix_constant_pH->return_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff,m_lambda_buff,N_buff);
     v_lambda_buff += dtf * a_lambda_buff;
     fix_constant_pH->reset_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff, m_lambda_buff);
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions
-----------------------------------------------------------------------*/

void FixNHConstantPH::nve_x()
{
  FixNH::nve_x();
  bigint ntimestep = update->ntimestep;
  fix_constant_pH->return_nparams(n_lambdas);
  fix_constant_pH->return_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);
  for (int i = 0; i < n_lambdas; i++)
     x_lambdas[i] += dtv * v_lambdas[i];
  if (lambda_integration_flags & BUFFER) {
     fix_constant_pH->return_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff,m_lambda_buff,N_buff);
     x_lambda_buff += dtv * v_lambda_buff;
     if (lambda_integration_flags & CONSTRAIN) constrain_lambdas();
     fix_constant_pH->reset_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff, m_lambda_buff);
  }
   
  fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);

  // This function sets the charges (qs) in the system based on the current value of x_lambdas and x_lambda_buffs
  fix_constant_pH->reset_qs();
}

/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of velocities
   ---------------------------------------------------------------------- */

void FixNHConstantPH::nh_v_temp()
{
  FixNH::nh_v_temp();
  // The timestep, the current step and the kT of course! 
  double dt = update->dt;
  bigint ntimestep = update->ntimestep;
  double kT = force->boltz * t_target;

  // Lets extract the parameters from the fix_constant_pH again
  fix_constant_pH->return_nparams(n_lambdas);
  fix_constant_pH->return_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);

  // and the buffer if present
  if (lambda_integration_flags & BUFFER)
     fix_constant_pH->return_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff,m_lambda_buff,N_buff);

  // Temperature
  double t_lambda_current;
  double t_lambda_target = t_target;
  fix_constant_pH->return_T_lambda(t_lambda_current);
  
  if (lambda_thermostat_type == LAMBDA_ANDERSEN) {
    double P = 1 - std::exp(-dt/t_andersen);
   

    if (which == NOBIAS) {
      // Dealing with lambdas
      for (int i = 0; i < n_lambdas; i++) {
        double r = static_cast<double>(rand()) / RAND_MAX;
        if (r < P) {
           double mean = 0.0;
           double sigma = std::sqrt(0.0019872041*4184.0*kT/ (10.0* m_lambdas[i]))/1000.0;
           v_lambdas[i] = random_normal(mean, sigma);
        }
        if (x_lambdas[i] < -0.1 || x_lambdas[i] > 1.1)
           v_lambdas[i] = -(x_lambdas[i]/std::abs(x_lambdas[i]))*std::abs(v_lambdas[i]);
      }
      // Dealing with the buffer
      if (lambda_integration_flags & BUFFER) {
         double r = static_cast<double>(rand())/ RAND_MAX;
         if (r < P) {
            double mean = 0.0;
            double sigma = std::sqrt(0.0019872041*4184.0*kT/ (10.0* m_lambda_buff))/1000.0;
            v_lambda_buff = random_normal(mean,sigma);
         }
         if (x_lambda_buff < -0.1 || x_lambda_buff > 1.1)
            v_lambda_buff = -(x_lambda_buff/std::abs(x_lambda_buff))*std::abs(v_lambda_buff);
      }
    } else if (which == BIAS) {
      // This needs to be implemented
      error->one(FLERR,"The bias keyword for the fix_nh_constant_pH has not been implemented yet!");
    }
  } else if (lambda_thermostat_type == LAMBDA_BOSSI) {
     //tau_t_bussi should be 1000
     
     double scaling_factor = std::sqrt(t_lambda_target/t_lambda_current);

     if (which == NOBIAS) {
        double friction = (t_lambda_current/t_lambda_target - 1.0) / tau_t_bussi;
        zeta_bussi += friction * dt; 

        // first, the lambdas
        for (int i = 0; i < n_lambdas; i++) {
           v_lambdas[i] *= scaling_factor;
           //v_lambdas[i] *= exp(-zeta*dt);
           if (x_lambdas[i] < -0.1 || x_lambdas[i] > 1.1)
              v_lambdas[i] = -(x_lambdas[i]/std::abs(x_lambdas[i]))*std::abs(v_lambdas[i]);
        }
        // and then the buffer 
        if (lambda_integration_flags & BUFFER) {
           v_lambda_buff *= scaling_factor;
           // v_lambda_buff *= exp(-zeta*dt);
           if (x_lambda_buff < -0.1 || x_lambda_buff > 1.1)
              v_lambda_buff = -(x_lambda_buff/std::abs(x_lambda_buff))*std::abs(v_lambda_buff);
        }
     } else if (which == BIAS) {
        // This needs to be implemented
        error->one(FLERR,"The bias keyword for the fix_nh_constant_pH has not been implemented yet!");
     }
  } else if (lambda_thermostat_type == LAMBDA_NOSEHOOVER) {  
     zeta_nose_hoover += dt * (t_lambda_current - t_lambda_target);

     if (which == NOBIAS) {
        // first the lambdas
        for (int i = 0; i < n_lambdas; i++) {
           v_lambdas[i] *= std::exp(-zeta_nose_hoover * dt);
           if (x_lambdas[i] < -0.1 || x_lambdas[i] > 1.1)
              v_lambdas[i] = -(x_lambdas[i]/std::abs(x_lambdas[i]))*std::abs(v_lambdas[i]);
        }
        // and then the buffer
        if (lambda_integration_flags & BUFFER) {
           v_lambda_buff *= std::exp(-zeta_nose_hoover * dt);
           if (x_lambda_buff < -0.1 || x_lambda_buff > 1.1)
              v_lambda_buff = -(x_lambda_buff/std::abs(x_lambda_buff))*std::abs(v_lambda_buff);
        }
     } else if (which == BIAS) {
        // This needs to be implemented
        error->one(FLERR,"The bias keyword for the fix_nh_constant_pH has not been implemented yet!");
     }
  }
  
  
  fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);

  if (lambda_integration_flags & BUFFER) fix_constant_pH->reset_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff, m_lambda_buff);
}

/* ---------------------------------------------------------------------
   applies the shake algorithm to the sum of the lambdas 

   It adds a  constraint according to the Donnine et al JCTC 2016 
   equation (13).

   The contraint equations were taken from the Tuckerman statistical mechanics
   book 2nd edition pages 106.

   Just one step of the shake iteration is enough as the constraint is very simple.
   
   --------------------------------------------------------------------- */

void FixNHConstantPH::constrain_lambdas()
{
   double sigma_lambda = 0.0;
   double sigma_mass_inverse  = 0.0;
   double domega; 

   for (int i = 0; i < n_lambdas; i++) {
      sigma_lambda += x_lambdas[i];
      sigma_mass_inverse += (1.0/m_lambdas[i]);
   }


   domega = -(mols_charge_change*sigma_lambda+buff_charge_change*N_buff*x_lambda_buff-total_charge)\ 
      / (mols_charge_change*sigma_mass_inverse + (static_cast<double>(N_buff*N_buff)*buff_charge_change*buff_charge_change/m_lambda_buff));

   for (int i = 0; i < n_lambdas; i++)
      x_lambdas[i] += (domega * mols_charge_change /m_lambdas[i]);

   x_lambda_buff += static_cast<double>(N_buff) * buff_charge_change * domega / m_lambda_buff;

}

/* ----------------------------------------------------------------------
   checking the total system charge after applying the constraint on the total charge
   ---------------------------------------------------------------------- */

void FixNHConstantPH::compute_q_total()
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

/* ----------------------------------------------------------------------
   random number generator
   ---------------------------------------------------------------------- */

double FixNHConstantPH::random_normal(double mean, double stddev)
{
  static std::mt19937 generator(std::random_device{}());
  std::normal_distribution<double> distribution(mean, stddev);
  return distribution(generator);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixNHConstantPH::memory_usage()
{
  double bytes = 0.0;
  bytes += 4.0*n_lambdas*sizeof(double); // x_lambdas, v_lambdas, a_lambdas and m_lambdas
  if (irregular) bytes += irregular->memory_usage();
  return bytes;
}
