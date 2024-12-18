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
    FixNH(lmp, narg-3, arg), fix_constant_pH_id(nullptr),
    x_lambdas(nullptr), v_lambdas(nullptr), a_lambdas(nullptr), m_lambdas(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);

  restart_global = 1;
  dynamic_group_allow = 1;
  time_integrate = 1;
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 0;
  ecouple_flag = 1;

  // default values

  pcouple = NONE;
  drag = 0.0;
  allremap = 1;
  id_dilate = nullptr;
  mtchain = mpchain = 3;
  nc_tchain = nc_pchain = 1;
  mtk_flag = 1;
  deviatoric_flag = 0;
  nreset_h0 = 0;
  eta_mass_flag = 1;
  omega_mass_flag = 0;
  etap_mass_flag = 0;
  flipflag = 1;
  dipole_flag = 0;
  dlm_flag = 0;
  p_temp_flag = 0;

  tcomputeflag = 0;
  pcomputeflag = 0;

  // turn on tilt factor scaling, whenever applicable

  dimension = domain->dimension;

  scaleyz = scalexz = scalexy = 0;
  if (domain->yperiodic && domain->xy != 0.0) scalexy = 1;
  if (domain->zperiodic && dimension == 3) {
    if (domain->yz != 0.0) scaleyz = 1;
    if (domain->xz != 0.0) scalexz = 1;
  }

  // set fixed-point to default = center of cell

  fixedpoint[0] = 0.5*(domain->boxlo[0]+domain->boxhi[0]);
  fixedpoint[1] = 0.5*(domain->boxlo[1]+domain->boxhi[1]);
  fixedpoint[2] = 0.5*(domain->boxlo[2]+domain->boxhi[2]);

  // used by FixNVTSllod to preserve non-default value

  mtchain_default_flag = 1;

  tstat_flag = 0;
  double t_period = 0.0;

  double p_period[6];
  for (int i = 0; i < 6; i++) {
    p_start[i] = p_stop[i] = p_period[i] = p_target[i] = 0.0;
    p_flag[i] = 0;
  }

  // process keywords

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} temp", style), error);
      tstat_flag = 1;
      t_start = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      t_target = t_start;
      t_stop = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      t_period = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (t_start <= 0.0 || t_stop <= 0.0)
        error->all(FLERR, "Target temperature for fix {} cannot be 0.0", style);
      iarg += 4;

    } else if (strcmp(arg[iarg],"iso") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} iso", style), error);
      pcouple = XYZ;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = p_period[1] = p_period[2] =
        utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg],"aniso") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} aniso", style), error);
      pcouple = NONE;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = p_period[1] = p_period[2] =
        utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg],"tri") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tri", style), error);
      pcouple = NONE;
      scalexy = scalexz = scaleyz = 0;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = p_period[1] = p_period[2] =
        utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      p_start[3] = p_start[4] = p_start[5] = 0.0;
      p_stop[3] = p_stop[4] = p_stop[5] = 0.0;
      p_period[3] = p_period[4] = p_period[5] =
        utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[3] = p_flag[4] = p_flag[5] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
        p_start[3] = p_stop[3] = p_period[3] = 0.0;
        p_flag[3] = 0;
        p_start[4] = p_stop[4] = p_period[4] = 0.0;
        p_flag[4] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} x", style), error);
      p_start[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = 1;
      deviatoric_flag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} y", style), error);
      p_start[1] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[1] = 1;
      deviatoric_flag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} z", style), error);
      p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[2] = 1;
      deviatoric_flag = 1;
      iarg += 4;
      if (dimension == 2) error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);

    } else if (strcmp(arg[iarg],"yz") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} yz", style), error);
      p_start[3] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[3] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[3] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[3] = 1;
      deviatoric_flag = 1;
      scaleyz = 0;
      iarg += 4;
      if (dimension == 2) error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);
    } else if (strcmp(arg[iarg],"xz") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} xz", style), error);
      p_start[4] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[4] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[4] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[4] = 1;
      deviatoric_flag = 1;
      scalexz = 0;
      iarg += 4;
      if (dimension == 2) error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);
    } else if (strcmp(arg[iarg],"xy") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} xy", style), error);
      p_start[5] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[5] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[5] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[5] = 1;
      deviatoric_flag = 1;
      scalexy = 0;
      iarg += 4;

    } else if (strcmp(arg[iarg],"couple") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} couple", style), error);
      if (strcmp(arg[iarg+1],"xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg+1],"xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg+1],"yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg+1],"xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg+1],"none") == 0) pcouple = NONE;
      else error->all(FLERR,"Illegal fix {} couple option: {}", style, arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"drag") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} drag", style), error);
      drag = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (drag < 0.0) error->all(FLERR, "Invalid fix {} drag argument: {}", style, drag);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ptemp") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} ptemp", style), error);
      p_temp = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_temp_flag = 1;
      if (p_temp <= 0.0) error->all(FLERR, "Invalid fix {} ptemp argument: {}", style, p_temp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} dilate", style), error);
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else {
        allremap = 0;
        delete[] id_dilate;
        id_dilate = utils::strdup(arg[iarg+1]);
        int idilate = group->find(id_dilate);
        if (idilate == -1)
          error->all(FLERR,"Fix {} dilate group ID {} does not exist", style, id_dilate);
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"tchain") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tchain", style), error);
      mtchain = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      // used by FixNVTSllod to preserve non-default value
      mtchain_default_flag = 0;
      if (mtchain < 1) error->all(FLERR, "Invalid fix {} tchain argument: {}", style, mtchain);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pchain") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} pchain", style), error);
      mpchain = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (mpchain < 0) error->all(FLERR, "Invalid fix {} pchain argument: {}", style, mpchain);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mtk") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} mtk", style), error);
      mtk_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tloop") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tloop", style), error);
      nc_tchain = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nc_tchain < 0) error->all(FLERR, "Invalid fix {} tloop argument: {}", style, nc_tchain);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ploop") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} ploop", style), error);
      nc_pchain = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nc_pchain < 0) error->all(FLERR, "Invalid fix {} ploop argument: {}", style, nc_pchain);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nreset") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} nreset", style), error);
      nreset_h0 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nreset_h0 < 0) error->all(FLERR, "Invalid fix {} nreset argument: {}", style, nreset_h0);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scalexy") == 0) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} scalexy", style), error);
      scalexy = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scalexz") == 0) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} scalexz", style), error);
      scalexz = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scaleyz") == 0) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} scaleyz", style), error);
      scaleyz = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"flip") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} flip", style), error);
      flipflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} update", style), error);
      if (strcmp(arg[iarg+1],"dipole") == 0) dipole_flag = 1;
      else if (strcmp(arg[iarg+1],"dipole/dlm") == 0) {
        dipole_flag = 1;
        dlm_flag = 1;
      } else error->all(FLERR, "Invalid fix {} update argument: {}", style, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"fixedpoint") == 0) {
      if (iarg+4 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} fixedpoint", style), error);
      fixedpoint[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      fixedpoint[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      fixedpoint[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    // disc keyword is also parsed in fix/nh/sphere

    } else if (strcmp(arg[iarg],"disc") == 0) {
      iarg++;

    // keywords erate, strain, and ext are also parsed in fix/nh/uef

    } else if (strcmp(arg[iarg],"erate") == 0) {
      iarg += 3;
    } else if (strcmp(arg[iarg],"strain") == 0) {
      iarg += 3;
    } else if (strcmp(arg[iarg],"ext") == 0) {
      iarg += 2;

    // keyword psllod is parsed in fix/nvt/sllod

    } else if (strcmp(arg[iarg],"psllod") == 0) {
      iarg += 2;

    } else if (strcmp(arg[iarg],"fix_constant_pH_id") == 0) {
       lambda_thermostat_type = NONE_LAMBDA;
       fix_constant_pH_id = utils::strdup(arg[iarg+1]);
       if (strcmp(arg[iarg+2],"none") == 0) {
          iarg+=2;
       }
       else if (strcmp(arg[iarg+2],"andersen") == 0) {
          lambda_thermostat_type = LAMBDA_ANDERSEN;
          t_andersen = utils::numeric(FLERR,arg[iarg+3],false,lmp);
          iarg+=3;
       }
       else if (strcmp(arg[iarg+2],"bossi") == 0) {
          lambda_thermostat_type = LAMBDA_BOSSI;
          iarg+=2;
       }
       else if (strcmp(arg[iarg+2],"nose-hoover") == 0) {
          lambda_thermostat_type = LAMBDA_NOSEHOOVER;
          Q_lambda_nose_hoover = utils::numeric(FLERR,arg[iarg+3],false,lmp);
          iarg+=3;
       }
       if (strcmp(arg[iarg],"buffer") == 0) {
          lambda_integration_flags |= BUFFER;
          iarg += 1;
       }
       if (strcmp(arg[iarg],"constrain_total_charge") == 0) {
          if (!(lambda_integration_flags & BUFFER))
             error->one(FLERR,"Constrain total charge in absence of a buffer is not supported yet!");
          lambda_integration_flags |= CONSTRAIN;
          mols_charge_change = utils::numeric(FLERR,arg[iarg+1],false,lmp);
          buff_charge_change = utils::numeric(FLERR,arg[iarg+2],false,lmp);
          total_charge = utils::numeric(FLERR,arg[iarg+3],false,lmp);
          iarg += 4;
       }

       // Here I cannot check if the constant_pH command has a buffer set as right now I have not set the fix_constant_pH pointer yet!!
    } else error->all(FLERR,"Unknown fix {} keyword: {}", style, arg[iarg]);
  }

  if (fix_constant_pH_id == nullptr) error->all(FLERR, "Invalid fix_nh constant_pH");

  // error checks

  if (dimension == 2 && (p_flag[2] || p_flag[3] || p_flag[4]))
    error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);
  if (dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);
  if (dimension == 2 && (scalexz == 1 || scaleyz == 1 ))
    error->all(FLERR,"Invalid fix {} command for a 2d simulation", style);

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == XYZ && dimension == 3 && p_flag[2] == 0)
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix {} command pressure settings", style);

  // require periodicity in tensile dimension

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR,"Cannot use fix {} on a non-periodic x dimension", style);
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR,"Cannot use fix {} on a non-periodic y dimension", style);
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR,"Cannot use fix {} on a non-periodic z dimension", style);

  // require periodicity in 2nd dim of off-diagonal tilt component

  if (p_flag[3] && domain->zperiodic == 0)
    error->all(FLERR, "Cannot use fix {} on a 2nd non-periodic dimension", style);
  if (p_flag[4] && domain->zperiodic == 0)
    error->all(FLERR, "Cannot use fix {} on a 2nd non-periodic dimension", style);
  if (p_flag[5] && domain->yperiodic == 0)
    error->all(FLERR, "Cannot use fix {} on a 2nd non-periodic dimension", style);

  if (scaleyz == 1 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot use fix {} with yz scaling when z is non-periodic dimension", style);
  if (scalexz == 1 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot use fix {} with xz scaling when z is non-periodic dimension", style);
  if (scalexy == 1 && domain->yperiodic == 0)
    error->all(FLERR,"Cannot use fix {} with xy scaling when y is non-periodic dimension", style);

  if (p_flag[3] && scaleyz == 1)
    error->all(FLERR,"Cannot use fix {} with both yz dynamics and yz scaling", style);
  if (p_flag[4] && scalexz == 1)
    error->all(FLERR,"Cannot use fix {} with both xz dynamics and xz scaling", style);
  if (p_flag[5] && scalexy == 1)
    error->all(FLERR,"Cannot use fix {} with both xy dynamics and xy scaling", style);

  if (!domain->triclinic && (p_flag[3] || p_flag[4] || p_flag[5]))
    error->all(FLERR,"Can not specify Pxy/Pxz/Pyz in fix {} with non-triclinic box", style);

  if (pcouple == XYZ && dimension == 3 &&
      (p_start[0] != p_start[1] || p_start[0] != p_start[2] ||
       p_stop[0] != p_stop[1] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[1] || p_period[0] != p_period[2]))
    error->all(FLERR,"Invalid fix {} pressure settings", style);
  if (pcouple == XYZ && dimension == 2 &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR,"Invalid fix {} pressure settings", style);
  if (pcouple == XY &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR,"Invalid fix {} pressure settings", style);
  if (pcouple == YZ &&
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2] ||
       p_period[1] != p_period[2]))
    error->all(FLERR,"Invalid fix {} pressure settings", style);
  if (pcouple == XZ &&
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[2]))
    error->all(FLERR,"Invalid fix {} pressure settings", style);

  if (dipole_flag) {
    if (strstr(style, "/sphere")) {
      if (!atom->omega_flag)
        error->all(FLERR,"Using update dipole flag requires atom attribute omega");
      if (!atom->radius_flag)
        error->all(FLERR,"Using update dipole flag requires atom attribute radius");
      if (!atom->mu_flag)
        error->all(FLERR,"Using update dipole flag requires atom attribute mu");
    } else {
      error->all(FLERR, "Must use a '/sphere' Nose-Hoover fix style for updating dipoles");
    }
  }

  if ((tstat_flag && t_period <= 0.0) ||
      (p_flag[0] && p_period[0] <= 0.0) ||
      (p_flag[1] && p_period[1] <= 0.0) ||
      (p_flag[2] && p_period[2] <= 0.0) ||
      (p_flag[3] && p_period[3] <= 0.0) ||
      (p_flag[4] && p_period[4] <= 0.0) ||
      (p_flag[5] && p_period[5] <= 0.0))
    error->all(FLERR,"Fix {} damping parameters must be > 0.0", style);

  // check that ptemp is not defined with a thermostat
  if (tstat_flag && p_temp_flag)
    error->all(FLERR,"Thermostat in fix {} is incompatible with ptemp command", style);

  // set pstat_flag and box change and restart_pbc variables

  pre_exchange_flag = 0;
  pstat_flag = 0;
  pstyle = ISO;

  for (int i = 0; i < 6; i++)
    if (p_flag[i]) pstat_flag = 1;

  if (pstat_flag) {
    if (p_flag[0]) box_change |= BOX_CHANGE_X;
    if (p_flag[1]) box_change |= BOX_CHANGE_Y;
    if (p_flag[2]) box_change |= BOX_CHANGE_Z;
    if (p_flag[3]) box_change |= BOX_CHANGE_YZ;
    if (p_flag[4]) box_change |= BOX_CHANGE_XZ;
    if (p_flag[5]) box_change |= BOX_CHANGE_XY;
    no_change_box = 1;
    if (allremap == 0) restart_pbc = 1;

    // pstyle = TRICLINIC if any off-diagonal term is controlled -> 6 dof
    // else pstyle = ISO if XYZ coupling or XY coupling in 2d -> 1 dof
    // else pstyle = ANISO -> 3 dof

    if (p_flag[3] || p_flag[4] || p_flag[5]) pstyle = TRICLINIC;
    else if (pcouple == XYZ || (dimension == 2 && pcouple == XY)) pstyle = ISO;
    else pstyle = ANISO;

    // pre_exchange only required if flips can occur due to shape changes

    if (flipflag && (p_flag[3] || p_flag[4] || p_flag[5]))
      pre_exchange_flag = pre_exchange_migrate = 1;
    if (flipflag && (domain->yz != 0.0 || domain->xz != 0.0 ||
                     domain->xy != 0.0))
      pre_exchange_flag = pre_exchange_migrate = 1;
  }

  // convert input periods to frequencies

  t_freq = 0.0;
  p_freq[0] = p_freq[1] = p_freq[2] = p_freq[3] = p_freq[4] = p_freq[5] = 0.0;

  if (tstat_flag) t_freq = 1.0 / t_period;
  if (p_flag[0]) p_freq[0] = 1.0 / p_period[0];
  if (p_flag[1]) p_freq[1] = 1.0 / p_period[1];
  if (p_flag[2]) p_freq[2] = 1.0 / p_period[2];
  if (p_flag[3]) p_freq[3] = 1.0 / p_period[3];
  if (p_flag[4]) p_freq[4] = 1.0 / p_period[4];
  if (p_flag[5]) p_freq[5] = 1.0 / p_period[5];

  // Nose/Hoover temp and pressure init

  size_vector = 0;

  if (tstat_flag) {
    int ich;
    eta = new double[mtchain];

    // add one extra dummy thermostat, set to zero

    eta_dot = new double[mtchain+1];
    eta_dot[mtchain] = 0.0;
    eta_dotdot = new double[mtchain];
    for (ich = 0; ich < mtchain; ich++) {
      eta[ich] = eta_dot[ich] = eta_dotdot[ich] = 0.0;
    }
    eta_mass = new double[mtchain];
    size_vector += 2*2*mtchain;
  }

  if (pstat_flag) {
    omega[0] = omega[1] = omega[2] = 0.0;
    omega_dot[0] = omega_dot[1] = omega_dot[2] = 0.0;
    omega_mass[0] = omega_mass[1] = omega_mass[2] = 0.0;
    omega[3] = omega[4] = omega[5] = 0.0;
    omega_dot[3] = omega_dot[4] = omega_dot[5] = 0.0;
    omega_mass[3] = omega_mass[4] = omega_mass[5] = 0.0;
    if (pstyle == ISO) size_vector += 2*2*1;
    else if (pstyle == ANISO) size_vector += 2*2*3;
    else if (pstyle == TRICLINIC) size_vector += 2*2*6;

    if (mpchain) {
      int ich;
      etap = new double[mpchain];

      // add one extra dummy thermostat, set to zero

      etap_dot = new double[mpchain+1];
      etap_dot[mpchain] = 0.0;
      etap_dotdot = new double[mpchain];
      for (ich = 0; ich < mpchain; ich++) {
        etap[ich] = etap_dot[ich] =
          etap_dotdot[ich] = 0.0;
      }
      etap_mass = new double[mpchain];
      size_vector += 2*2*mpchain;
    }

    if (deviatoric_flag) size_vector += 1;
  }

  if (pre_exchange_flag) irregular = new Irregular(lmp);
  else irregular = nullptr;

  // initialize vol0,t0 to zero to signal uninitialized
  // values then assigned in init(), if necessary

  vol0 = t0 = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNHConstantPH::~FixNHConstantPH()
{
  FixNH::~FixNH();
  if (fix_constant_pH_id) delete [] fix_constant_pH_id;
  if (x_lambdas) delete [] x_lambdas;
  if (v_lambdas) delete [] v_lambdas;
  if (a_lambdas) delete [] a_lambdas;
  if (m_lambdas) delete [] m_lambdas; 
}

/* ---------------------------------------------------------------------- */

void FixNHConstantPH::init()
{
  FixNH::init();
  fix_constant_pH = static_cast<FixConstantPH*>(modify->get_fix_by_id(fix_constant_pH_id));
  fix_constant_pH->return_nparams(n_lambdas);
   
  x_lambdas = new double[n_lambdas];
  v_lambdas = new double[n_lambdas];
  a_lambdas = new double[n_lambdas];
  m_lambdas = new double[n_lambdas];

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
     if (lambda_integration_flags & CONSTRAIN) contrain_lambdas();
     fix_constant_pH->reset_buff_params(x_lambda_buff,v_lambda_buff,a_lambda_buff, m_lambda_buff);
  }
   
  fix_constant_pH->reset_params(x_lambdas,v_lambdas,a_lambdas,m_lambdas);
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

void FixNHConstantPH::contrain_lambdas()
{
   double sigma_lambda = 0.0;
   double sigma_mass_inverse  = 0.0;
   double domega; 


   for (int i = 0; i < n_lambdas; i++) {
      sigma_lambda += x_lambdas[i];
      sigma_mass_inverse += (1.0/m_lambdas[i]);
   }

   mols_charge_change, buff_charge_change, total_charge;
   domega = -(mols_charge_change*sigma_lambda+buff_charge_change*N_buff*x_lambda_buff-total_charge)\ 
      / (mols_charge_change*sigma_mass_inverse + (static_cast<double>(N_buff*N_buff)*buff_charge_change*buff_charge_change/m_lambda_buff));

   for (int i = 0; i < n_lambdas; i++)
      x_lambdas[i] += (domega * mols_charge_change /m_lambdas[i]);

   x_lambda_buff += static_cast<double>(N_buff) * buff_charge_change * domega / m_lambda_buff;

}

/* ----------------------------------------------------------------------
   random number generator
   -----------------------------------------------------------------------*/

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
