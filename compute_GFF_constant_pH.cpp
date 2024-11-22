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
       protonable(nullptr), typePerProtMol(nullptr), pH1qs(nullptr), pH2qs(nullptr)
{
   if (narg < 6) error->all(FLERR, "Illegal number of arguments in compute GFF constant pH");
   // Reading the file that contains the charges before and after protonation/deprotonation
   if (comm->me == 0) {
      pHStructureFile = fopen(arg[3],"r"); // The command reads the file the type and charge of each atom before and after protonation
      if (pHStructureFile == nullptr)
         error->all(FLERR,"Unable to open the file");
   }
   lambda = utils::numeric(FLERR,arg[4],false,lmp);
   dlambda = utils::numeric(FLERR,arg[5],false,lmp);
   if (dlambda < 0) error->all(FLERR,"Illegal compute GFF constant pH dlambda value {}", dlambda);
   typeHW = utils::inumeric(FLERR,arg[6],false,lmp);
   if (typeHW > atom->ntypes) error->all(FLERR,"Illegal compute GFF constant pH HW atom type {}", typeHW);

   peflag = 1;
   peatomflag = 1;
   peratom_flag = 1;
   comm_reverse = 1;
    
    
   scalar_flag = 0;
   vector_flag = 1;
   size_vector = 4;
   peratom_flag = 1; // I need to have per atom energies tallied. 
    
   extvector = 0;

   vector = new double[size_vector];

   /* Since I have a deallocate_storage in the the destructor
      I need this one so that if the number of steps is zero 
      there is some allocated memory to get deallocated.
   */
   allocate_storage();
}

/* ---------------------------------------------------------------------------- */
ComputeGFFConstantPH::~ComputeGFFConstantPH()
{
   delete [] vector;
   deallocate_storage();
   if (protonable) memory->destroy(protonable);
   if (typePerProtMol) memory->destroy(typePerProtMol);
   if (pH1qs) memory->destroy(pH1qs);
   if (pH2qs) memory->destroy(pH2qs);
}

/* ---------------------------------------------------------------------------- */
void ComputeGFFConstantPH::init()
{
   // Reading the structure of protonable states before and after protonation.
   read_pH_structure_files();
}

/* ---------------------------------------------------------------------------- */
void ComputeGFFConstantPH::setup()
{
   // Checking if we have enough hydronium ions to neutralize the system
   check_num_HWs();
   // Calculating the dq
   calculate_dq();
}

/* --------------------------------------------------------------------- */
void ComputeGFFConstantPH::compute_vector()
{
   // Calculating the HA, HB, HC and dH_dLambda
   compute_Hs();
   
   vector[0] = HA;
   vector[1] = HB;
   vector[2] = HC;
   vector[3] = dH_dLambda;
}

/* ---------------------------------------------------------------------- */

void ComputeGFFConstantPH::read_pH_structure_files()
{
   /* File format
    * Comment 
    * pHnTypes
    * type1,  number of type1 atoms in the protonable molecule, qBeforeProtonation, qAfterProtonation
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

/* ---------------------------------------------------------------------- */

void ComputeGFFConstantPH::calculate_dq()
{
   double q_total_1 = 0.0;
   double q_total_2 = 0.0;

   int ntypes = atom->ntypes;

   for (int i = 1; i < ntypes+1; i++)
   {
       q_total_1 += typePerProtMol[i] * protonable[i] * pH1qs[i]; // if it is not protonable the protonable[i] == 0
       q_total_2 += typePerProtMol[i] * protonable[i] * pH2qs[i];
   }
   
   dq = (q_total_2 - q_total_1);
}

/* ----------------------------------------------------------------------
   Checking if we have enough of HWs for neutralizing the system total charge
   ---------------------------------------------------------------------- */

void ComputeGFFConstantPH::check_num_HWs()
{
   double tol = 1e-5;
   int * type = atom->type;
   int nlocal = atom->nlocal;
   int num_local, num;
   num = 0;
   num_local = 0;
   for (int i = 0; i < nlocal; i++)
   {
      if (type[i] == typeHW)
        num_local++;
   }

   MPI_Allreduce(&num_local,&num,1,MPI_INT,MPI_SUM,world);
   num_HWs = num;
   if (comm->me == 0) error->warning(FLERR,"The num_HWs = {}",num_HWs);
}

/* ----------------------------------------------------------------------- */

void ComputeGFFConstantPH::compute_Hs()
{
  if (nmax > atom->nmax) {
     nmax = atom->nmax;
     deallocate_storage();
     allocate_storage();
  }
  HC = compute_epair();
  backup_restore_qfev<1>();      // backup charge, force, energy, virial array values
  modify_lambda(lambda-dlambda); //should define a change_parameters(const int);
  update_lmp(); // update the lammps force and virial values
  HA = compute_epair(); 
  backup_restore_qfev<-1>();        // restore charge, force, energy, virial array values
  modify_lambda(lambda+dlambda); //should define a change_parameters(const double);
  update_lmp();
  HB = compute_epair();           // HB is for the protonated state with lambda==1 
  backup_restore_qfev<-1>();      // restore charge, force, energy, virial array values
  dH_dLambda = (HB-HA)/(2*dlambda);
}

/* ----------------------------------------------------------------------
   manage storage for charge, force, energy, virial arrays
   taken from src/FEP/compute_fep.cpp
------------------------------------------------------------------------- */

void ComputeGFFConstantPH::allocate_storage()
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

void ComputeGFFConstantPH::deallocate_storage()
{
  memory->destroy(q_orig);
  memory->destroy(f_orig);
  memory->destroy(peatom_orig);
  memory->destroy(pvatom_orig);
  if (force->kspace)
  {
      memory->destroy(keatom_orig);
      memory->destroy(kvatom_orig);
  }

  f_orig = nullptr;
  peatom_orig = keatom_orig = nullptr;
  pvatom_orig = kvatom_orig = nullptr;
}

/* ----------------------------------------------------------------------
   Forward-reverse copy function to be used in backup_restore_qfev()
   ---------------------------------------------------------------------- */

template  <int direction>
void ComputeGFFConstantPH::forward_reverse_copy(double& a, double& b)
{
    if (direction == 1) a = b;
    if (direction == -1) b = a;
}

template  <int direction>
void ComputeGFFConstantPH::forward_reverse_copy(double* a, double* b, int m)
{
    for (int i = 0; i < m; i++) forward_reverse_copy<direction>(a[i],b[i]);
}

template  <int direction>
void ComputeGFFConstantPH::forward_reverse_copy(double** a, double** b, int m, int n)
{
    for (int i = 0; i < m; i++) forward_reverse_copy<direction>(a[i],b[i],n);
}

/* ----------------------------------------------------------------------
   backup and restore arrays with charge, force, energy, virial
   taken from src/FEP/compute_fep.cpp
   backup ==> direction == 1
   restore ==> direction == -1
------------------------------------------------------------------------- */

template <int direction>
void ComputeGFFConstantPH::backup_restore_qfev()
{
    int i;

    int nall = atom->nlocal + atom->nghost;
    int natom = atom->nlocal;
    if (force->newton || force->kspace->tip4pflag) natom += atom->nghost;

    double** f = atom->f;
    forward_reverse_copy<direction>(f_orig, f, natom, 3);

    double* q = atom->q;
    forward_reverse_copy<direction>(q_orig, q, natom);

    forward_reverse_copy<direction>(eng_vdwl_orig, force->pair->eng_vdwl);
    forward_reverse_copy<direction>(eng_coul_orig, force->pair->eng_coul);

    
    forward_reverse_copy<direction>(pvirial_orig, force->pair->virial, 6);

    if (update->eflag_atom) {
        double* peatom = force->pair->eatom;
        forward_reverse_copy<direction>(peatom_orig, peatom, natom);
    }
    if (update->vflag_atom) {
        double** pvatom = force->pair->vatom;
        forward_reverse_copy<direction>(pvatom_orig, pvatom, natom, 6);
    }

    if (force->kspace) {
        forward_reverse_copy<direction>(energy_orig, force->kspace->energy);
        forward_reverse_copy<direction>(kvirial_orig, force->kspace->virial, 6);

        if (update->eflag_atom) {
            double* keatom = force->kspace->eatom;
            forward_reverse_copy<direction>(keatom_orig, keatom, natom);
        }
        if (update->vflag_atom) {
            double** kvatom = force->kspace->vatom;
            forward_reverse_copy<direction>(kvatom_orig, kvatom, natom, 6);
        }
    }
}

/* --------------------------------------------------------------

   -------------------------------------------------------------- */
   
void ComputeGFFConstantPH::modify_lambda(const double& scale)
{
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    int * type = atom->type;
    int ntypes = atom->ntypes;
    double * q = atom->q;


    double q_changes_local[4] = {0.0,0.0,0.0,0.0};
    double q_changes[4] = {0.0,0.0,0.0,0.0};

	
    // update the charges
    for (int i = 0; i < nlocal; i++)
    {
       if (protonable[type[i]] == 1)
       {
           double q_init = q[i];
           q[i] = pH1qs[type[i]] + scale * (pH2qs[type[i]] - pH1qs[type[i]]); // scale == 1 should be for the protonated state
	   q_changes_local[0]++;
	   q_changes_local[2] += (q[i] - q_init);
       }
       if (type[i] == typeHW)
       {
	   double q_init = q[i];
	   q[i] = q_init + (-scale) * dq / static_cast<double> (num_HWs); //The total charge should be neutral
	   q_changes_local[1]++;
	   q_changes_local[3] += (q[i] - q_init);
       }
    }

    /* The purpose of this part this is just to debug the total charge.
       So, in the final version of the code this part should be 
       commented out!
    */
    /*MPI_Allreduce(&q_changes_local,&q_changes,4,MPI_DOUBLE,MPI_SUM,world);
    if (comm->me == 0) error->warning(FLERR,"protonable q change = {}, HW q change = {}, protonable charge change = {}, HW charge change = {}",q_changes[0],q_changes[1],q_changes[2],q_changes[3]);
    */
    compute_q_total();
}

/* ----------------------------------------------------------------------
   modify force and kspace in lammps according
   ---------------------------------------------------------------------- */

void ComputeGFFConstantPH::update_lmp() {
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
    //if (fixgpu) fixgpu->post_force(vflag);
}

/* --------------------------------------------------------------------- */

void ComputeGFFConstantPH::compute_q_total()
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

double ComputeGFFConstantPH::compute_epair()
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
