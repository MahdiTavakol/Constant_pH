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


#include "compute_GFF_constant_pH.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "timer.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------------- */

ComputeGFFConstantPH::ComputeGFFConstantPH(LAMMPS* lmp, int narg, char** arg) : Compute(lmp, narg, arg), 
       fix_constant_pH_id(nullptr)
{
   if (narg < 4) error->all(FLERR, "Illegal number of arguments in compute GFF constant pH");
   lambda = utils::numeric(FLERR,arg[3],false,lmp);
   dlambda = utils::numeric(FLERR,arg[4],false,lmp);
   if (dlambda < 0) error->all(FLERR,"Illegal compute GFF constant pH dlambda value {}", dlambda);

    
    
   array_flag = 1;
   size_array_cols = 4;
   size_array_rows_variable = 1;
    
   extvector = 0;

}

/* ---------------------------------------------------------------------------- */
ComputeGFFConstantPH::~ComputeGFFConstantPH()
{
   if (fix_constant_pH_id) delete [] fix_constant_pH_id;
   deallocate_storage();
}

/* ---------------------------------------------------------------------------- */
void ComputeGFFConstantPH::init()
{
  fix_constant_pH = static_cast<FixConstantPH*>(modify->get_fix_by_id(fix_constant_pH_id));
  fix_constant_pH->return_nparams(n_lambdas);

  allocated_storage();
}

/* ---------------------------------------------------------------------------- */
void ComputeGFFConstantPH::setup()
{
   
}

/* --------------------------------------------------------------------- */
void ComputeGFFConstantPH::compute_array()
{
   // Checking if there is enough space in our arrays for all the lambas
   double n_params;
   fix_constant_pH->return_nparams(n_params);
   if (n_lambdas < n_params) {
      n_lambdas = n_params;
      deallocate_storage();
      allocate_storage();
   }
	
   // Calculating the HA, HB, HC and dH_dLambda
   for (int i = 0; i < n_lambda; i++) {
      fix_constant_pH->return_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);
      std::fill(x_lambdas,x_lambdas+n_lambdas,lambda);

      
      fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);
      fix_constant_pH->calculate_H_once();
      fix_constant_pH->return_H_lambdas(HCs);
	   
      x_lambdas[i] = lambda - dlambda;
      fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);
      fix_constant_pH->calculate_H_once();
      fix_constant_pH->return_H_lambas(HAs);

      x_lambdas[i] = lambda + dlambda;
      fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);
      fix_constant_pH->calculate_H_once();
      fix_constant_pH->return_H_lambdas(HBs);


      for (int j = 0; j < n_lambdas; j++) {
         array[j][0] = HAs[j];
	 array[j][1] = HBs[j];
	 array[j][2] = HCs[j];
	 array[j][0] = (HBs[j]-HAs[j])/(2*dlambda);
      }
   }
}

/* ----------------------------------------------------------------------
   manage storage for charge, force, energy, virial arrays
   taken from src/FEP/compute_fep.cpp
------------------------------------------------------------------------- */

void ComputeGFFConstantPH::allocate_storage()
{
  memory->create(array,nlambdas,4,"compute_GFF_constant_pH:array");
	
  memory->create(x_lambdas, nlambdas, 1, "compute_GFF_constant_pH:x_lambdas");
  memory->create(v_lambdas, nlambdas, 1, "compute_GFF_constant_pH:v_lambdas");
  memory->create(a_lambdas, nlambdas, 1, "compute_GFF_constant_pH:a_lambdas");
  memory->create(m_lambdas, nlambdas, 1, "compute_GFF_constant_pH:m_lambdas");
  memory->create(H_lambdas, nlambdas, 1, "compute_GFF_constant_pH:H_lambdas");
	
  memory->create(HAs, nlambdas, 1, "compute_GFF_constant_pH:HAs");
  memory->create(HBs, nlambdas, 1, "compute_GFF_constant_pH:HBs");
  memory->create(HCs, nlambdas, 1, "compute_GFF_constant_pH:HCs");
}

/* ---------------------------------------------------------------------- */

void ComputeGFFConstantPH::deallocate_storage()
{
  if (array) memory->destroy(array);
	
  if (x_lambdas) memory->destory(x_lambdas);
  if (v_lambdas) memory->destory(v_lambdas);
  if (a_lambdas) memory->destory(a_lambdas);
  if (m_lambdas) memory->destory(m_lambdas);
  if (H_lambdas) memory->destory(H_lambdas);
	
  if (HAs) memory->destroy(HAs);
  if (HBs) memory->destroy(HBs);
  if (HCs) memory->destroy(HCs);

  array = nullptr;

  x_lambdas = nullptr;
  v_lambdas = nullptr;
  a_lambdas = nullptr;
  m_lambdas = nullptr;
  H_lambdas = nullptr;
	
  HAs = nullptr;
  HBs = nullptr;
  HCs = nullptr; 
}
