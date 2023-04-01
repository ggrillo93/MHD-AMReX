
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

std::function<double(double)> limFunc; // slope limiter for MUSCL-Hancock
int orient; // allows for different tests
int order = 2; // order = 2 uses MUSCL-Hancock, otherwise accuracy is first order

int AmrLevelAdv::verbose = 0;
Real AmrLevelAdv::cfl = 0.8; // Default value - can be overwritten in settings file
int AmrLevelAdv::do_reflux = 1;

State wR0, wL0; // initial states

State varLim;

int AmrLevelAdv::NUM_STATE = wR0.size();
int AmrLevelAdv::NUM_GROW = 2;  // number of ghost cells

// Default constructor.  Builds invalid object.
AmrLevelAdv::AmrLevelAdv()
{
  // Flux registers store fluxes at patch boundaries to ensure fluxes are conservative between AMR levels
  flux_reg = 0;
}

// The basic constructor.
AmrLevelAdv::AmrLevelAdv(Amr &papa, int lev, const Geometry &level_geom, const BoxArray &bl, const DistributionMapping &dm , Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time)
{
  // Flux registers are only required if AMR is actually used, and if flux fix up is being done (recommended)
  flux_reg = 0;
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
  }
}

// The destructor.
AmrLevelAdv::~AmrLevelAdv()
{
  delete flux_reg;
}

// Restart from a checkpoint file.
void AmrLevelAdv::restart(Amr &papa, std::istream &is, bool bReadSpecial)
{
  AmrLevel::restart(papa, is, bReadSpecial);

  BL_ASSERT(flux_reg == 0);
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
  }
}

// Write a checkpoint file - format is handled automatically by AMReX
void AmrLevelAdv::checkPoint(const std::string &dir, std::ostream &os, VisMF::How how, bool dump_old)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

// Write a plotfile to specified directory - format is handled automatically by AMReX.
void AmrLevelAdv::writePlotFile(const std::string &dir, std::ostream &os, VisMF::How how)
{
  AmrLevel::writePlotFile(dir, os, how);
}

// Define data descriptors.
void AmrLevelAdv::variableSetUp()
{
  BL_ASSERT(desc_lst.size() == 0);

  read_params();

  const int storedGhostZones = 0;

  // Setting up a container for a variable, or vector of variables:
  // Phi_Type: Enumerator for this variable type
  // IndexType::TheCellType(): AMReX can support cell-centred and vertex-centred variables (cell centred here)
  // StateDescriptor::Point: Data can be a point in time, or an interval over time (point here)
  // storedGhostZones: Ghost zones can be stored (e.g. for output).  Generally set to zero.
  // NUM_STATE: Number of variables in the variable vector
  // cell_cons_interp: Controls interpolation between levels - cons_interp is good for finite volume
  desc_lst.addDescriptor(Phi_Type, IndexType::TheCellType(), StateDescriptor::Point, storedGhostZones, NUM_STATE, &cell_cons_interp);

  // Set up boundary condition
  int lo_bc[amrex::SpaceDim];
  int hi_bc[amrex::SpaceDim];
  // AMReX has pre-set BCs, including periodic (int_dir) and transmissive (foextrap)
  for (int i = 0; i < amrex::SpaceDim; ++i)
  {
    lo_bc[i] = hi_bc[i] = BCType::foextrap; // need to make this dependent on the settings file
  }

  // Object for storing all the boundary conditions
  BCRec bc(lo_bc, hi_bc);

  // Set up variable-specific information; needs to be done for each variable in NUM_STATE
  // Phi_Type: Enumerator for the variable type being set
  // 0: Position of the variable in the variable vector.  Single variable for advection.
  // phi: Name of the variable - appears in output to identify what is being plotted
  // bc: Boundary condition object for this variable (defined above)
  // BndryFunc: Function for setting boundary conditions.  For basic BCs, AMReX can handle these automatically
  desc_lst.setComponent(Phi_Type, 0, "rho", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 1, "mom_x", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 2, "mom_y", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 3, "mom_z", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 4, "energy", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 5, "Bx", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 6, "By", bc, StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 7, "Bz", bc, StateDescriptor::BndryFunc(nullfill));
  // desc_lst.setComponent(Phi_Type, 8, "psi", bc, StateDescriptor::BndryFunc(nullfill)); // divergence cleaning variable
}

// Cleanup data descriptors at end of run.
void AmrLevelAdv::variableCleanUp()
{
  desc_lst.clear();
}

template<typename T>
void AmrLevelAdv::calculateFluxes(const int i, const int j, const int k, const int iOffset, const int jOffset, const int kOffset, const Dim3 lo, State& uLL, State& uL, State& u0, State& uR, State& uLRe_L, State& uLRe_R, State& u0Re_L, State& u0Re_R, State& fluxVec, const T& arr, Real step, const unsigned int d) {
  /* Calculates intercell flux at the left boundary of cell u0 */
  if ((d == 0 && i == lo.x) || (d == 1 && j == lo.y)) { 
    for (int ii = 0; ii < NUM_STATE; ii++) {
      uLL[ii] = arr(i - 2 * iOffset, j - 2 * jOffset, k - 2 * kOffset, ii);
      uL[ii] = arr(i - iOffset, j - jOffset, k - kOffset, ii);
      u0[ii] = arr(i, j, k, ii);
      uR[ii] = arr(i + iOffset, j + jOffset, k + kOffset, ii);
    }
    reconstruct(uLL, uL, u0, uLRe_L, uLRe_R, varLim);
    reconstruct(uL, u0, uR, u0Re_L, u0Re_R, varLim);
    if (order == 2) {
      halfTimeEvol(uLRe_L, uLRe_R, step, d);
      halfTimeEvol(u0Re_L, u0Re_R, step, d);
    }
  }
  else {
    uLRe_L = u0Re_L;
    uLRe_R = u0Re_R;
    uL = u0;
    u0 = uR;
    for (int ii = 0; ii < NUM_STATE; ii++) {
      uR[ii] = arr(i + iOffset, j + jOffset, k + kOffset, ii);
    }
    reconstruct(uL, u0, uR, u0Re_L, u0Re_R, varLim);
    if (order == 2) {
      halfTimeEvol(u0Re_L, u0Re_R, step, d);
    }
  }
  HLLCFlux(uLRe_R, u0Re_L, fluxVec, d);
}

// Initialize grid data at problem start-up.
void AmrLevelAdv::initData()
{
  
  // Loop over grids, call FORTRAN function to init with data.
  const Real *dx = geom.CellSize();
  // Position of the bottom left corner of the domain
  const Real *prob_lo = geom.ProbLo();
  // Create a multifab which can store the initial data
  MultiFab &S_new = get_new_data(Phi_Type);
  Real cur_time = state[Phi_Type].curTime();

  if (verbose)
  {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

  // Ensure uninitialised data is not used.
  const Real dX = dx[0];
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);
  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);

  const Real probLoX = prob_lo[0];
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);
  const Real probLoZ = (amrex::SpaceDim > 2 ? prob_lo[2] : 0.0);

  Real var;
  State uL0, uR0;
  conservative(wL0, uL0);
  conservative(wR0, uR0);

  // Loop over all the patches at this level
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto &arr = S_new.array(mfi);

    for (int k = lo.z; k <= hi.z; k++)
    {
      const Real z = probLoZ + (double(k) + 0.5) * dZ;
      for (int j = lo.y; j <= hi.y; j++)
      {
        const Real y = probLoY + (double(j) + 0.5) * dY;
        for (int i = lo.x; i <= hi.x; i++)
        {
          const Real x = probLoX + (double(i) + 0.5) * dX;
          if (orient == 1) { // divide states along line x = 0.5
            var = x;
          }
          else if (orient == 2) { // divide states along line y = 0.5
            var = y;
          }
          else if (orient == 3) { // divide states along line x + y = 0.5
            var = x + y - 0.5;
          }
          else if (orient == 4) { // divide states along circle of radius 1
            var = (x - 1) * (x - 1) + (y - 1) * (y - 1) + 0.34; // cylindrical explosion
          }
          if (var < 0.5) {
            for (int ii = 0; ii < NUM_STATE; ii++) {
              arr(i, j, k, ii) = uL0[ii];
            }
          }
          else {
            for (int ii = 0; ii < NUM_STATE; ii++) {
              arr(i, j, k, ii) = uR0[ii];
            }
          }
        }
      }
    }
  }

  if (verbose)
  {
    amrex::Print() << "Done initializing the level " << level << " data " << std::endl;
  }
}

// Initialize data on this level from another AmrLevelAdv (during regrid).
void AmrLevelAdv::init(AmrLevel &old)
{

  AmrLevelAdv *oldlev = (AmrLevelAdv *)&old;
  
  // Create new grid data by fillpatching from old.
  Real dt_new = parent->dtLevel(level);
  Real cur_time = oldlev->state[Phi_Type].curTime();
  Real prev_time = oldlev->state[Phi_Type].prevTime();
  Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  MultiFab &S_new = get_new_data(Phi_Type);

  const int zeroGhosts = 0;
  // FillPatch takes the data from the first argument (which contains
  // all patches at a refinement level) and fills (copies) the
  // appropriate data onto the patch specified by the second argument:
  // old: Source data
  // S_new: destination data
  // zeroGhosts: If this is non-zero, ghost zones could be filled too - not needed for init routines
  // cur_time: AMReX can attempt interpolation if a different time is specified - not recommended for advection eq.
  // Phi_Type: Specify the type of data being set
  // 0: This is the first data index that is to be copied
  // NUM_STATE: This is the number of states to be copied
  FillPatch(old, S_new, zeroGhosts, cur_time, Phi_Type, 0, NUM_STATE);
}

// Initialize data on this level after regridding if old level did not previously exist
void AmrLevelAdv::init()
{
  Real dt = parent->dtLevel(level);
  Real cur_time = getLevel(level - 1).state[Phi_Type].curTime();
  Real prev_time = getLevel(level - 1).state[Phi_Type].prevTime();

  Real dt_old = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level - 1);

  setTimeLevel(cur_time, dt_old, dt);
  MultiFab &S_new = get_new_data(Phi_Type);

  FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

// Advance grids at this level in time.
Real AmrLevelAdv::advance(Real time, Real dt, int iteration, int ncycle)
{

  MultiFab &S_mm = get_new_data(Phi_Type);

  for (int k = 0; k < NUM_STATE_TYPE; k++)
  {
    state[k].allocOldData();
    state[k].swapTimeLevels(dt);
  }

  // S_new is the MultiFab that will be operated upon to update the data
  MultiFab &S_new = get_new_data(Phi_Type);

  const Real prev_time = state[Phi_Type].prevTime();
  const Real cur_time = state[Phi_Type].curTime();
  const Real ctr_time = 0.5 * (prev_time + cur_time);

  const Real *dx = geom.CellSize();
  const Real *prob_lo = geom.ProbLo();

  // Get pointers to Flux registers, or set pointer to zero if not there.
  FluxRegister *fine = 0;
  FluxRegister *current = 0;

  int finest_level = parent->finestLevel();

  // If we are not on the finest level, fluxes may need correcting
  // from those from finer levels.  To start this process, we set the
  // flux register values to zero
  if (do_reflux && level < finest_level)
  {
    fine = &getFluxReg(level + 1);
    fine->setVal(0.0);
  }

  // If we are not on the coarsest level, the fluxes are going to be
  // used to correct those on coarser levels.  We get the appropriate
  // flux level to include our fluxes within
  if (do_reflux && level > 0)
  {
    current = &getFluxReg(level);
  }

  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];

  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  // storage is increased by 1 cell in the direction of the flux.
  // This is only needed if refluxing is happening, otherwise fluxes
  // don't need to be stored, just used
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_new.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  
  FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  // Fill periodic boundaries where they exist.  More accurately, the
  // FillBoundary call will fill overlapping boundaries (with periodic
  // domains effectively being overlapping).  It also takes care of
  // AMR patch and CPU boundaries.
  Sborder.FillBoundary(geom.periodicity());

  Vector<BCRec> bc(NUM_STATE); // copied from AmReX documentation
  for (int n = 0; n < NUM_STATE; ++n)
  {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      {
          if (geom.isPeriodic(idim))
          {
              bc[n].setLo(idim, BCType::int_dir); // interior
              bc[n].setHi(idim, BCType::int_dir);
          }
          else
          {
              bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation
              bc[n].setHi(idim, BCType::foextrap);
          }
      }
  }

  FillDomainBoundary(Sborder, geom, bc);

  // Initialize necessary state vectors
  State fluxVec, uLL, uL, u0, uR, uLRe_L, uLRe_R, u0Re_L, u0Re_R;

  for (unsigned int d = 0; d < amrex::SpaceDim; d++)
  {
    Real step = dt / dx[d];

    const int iOffset = (d == 0 ? 1 : 0);
    const int jOffset = (d == 1 ? 1 : 0);
    const int kOffset = (d == 2 ? 1 : 0);

    // Loop over all the patches at this level

    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box &bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto &arr = Sborder.array(mfi);
      const auto &fluxArr = fluxes[d].array(mfi);
    
      if (d == 0) {
        for (int k = lo.z; k <= hi.z + kOffset; k++) {
          for (int j = lo.y; j <= hi.y + jOffset; j++) {
            for (int i = lo.x; i <= hi.x + iOffset; i++) {
              calculateFluxes(i, j, k, iOffset, jOffset, kOffset, lo, uLL, uL, u0, uR, uLRe_L, uLRe_R, u0Re_L, u0Re_R, fluxVec, arr, step, d);
              for (int ii = 0; ii < NUM_STATE; ii++) {
                fluxArr(i, j, k, ii) = fluxVec[ii];
              }
            }
          }
        }
      }
      else { // need to reverse the loops to be able to copy cell states in our flux calculation. Will need to do this again for 3D
        for (int k = lo.z; k <= hi.z + kOffset; k++) {
          for (int i = lo.x; i <= hi.x + iOffset; i++) {
            for (int j = lo.y; j <= hi.y + jOffset; j++) {
              calculateFluxes(i, j, k, iOffset, jOffset, kOffset, lo, uLL, uL, u0, uR, uLRe_L, uLRe_R, u0Re_L, u0Re_R, fluxVec, arr, step, d);
              for (int ii = 0; ii < NUM_STATE; ii++) {
                fluxArr(i, j, k, ii) = fluxVec[ii];
              }
            }
          }
        }
      }
    
      for (int k = lo.z; k <= hi.z; k++) {
        for (int j = lo.y; j <= hi.y; j++) {
          for (int i = lo.x; i <= hi.x; i++) {
            for (int ii = 0; ii < NUM_STATE; ii++) {
              arr(i, j, k, ii) = arr(i, j, k, ii) - step * (fluxArr(i + iOffset, j + jOffset, k + kOffset, ii) - fluxArr(i, j, k, ii));
            }
          }
        }
      }
    }

    // We need to compute boundary conditions again after each update
    FillDomainBoundary(Sborder, geom, bc);

    Sborder.FillBoundary(geom.periodicity());

    // The fluxes now need scaling for the reflux command.
    // This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    if (do_reflux)
    {
      Real scaleFactor = dt;
      for (int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
        // Fluxes don't need scaling by dx[d]
        if (scaledir == d)
        {
          continue;
        }
        scaleFactor *= dx[scaledir];
      }
      // The mult function automatically multiplies entries in a multifab by a scalar
      // scaleFactor: The scalar to multiply by
      // 0: The first data index in the multifab to multiply
      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
    }
  }
  // The updated data is now copied to the S_new multifab.  This means
  // it is now accessible through the get_new_data command, and AMReX
  // can automatically interpolate or extrapolate between layers etc.
  // S_new: Destination
  // Sborder: Source
  // Third entry: Starting variable in the source array to be copied (the zeroth variable in this case)
  // Fourth entry: Starting variable in the destination array to receive the copy (again zeroth here)
  // NUM_STATE: Total number of variables being copied
  // Sixth entry: Number of ghost cells to be included in the copy (zero in this case, since only real data is needed for S_new)
  MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, 0);

  // Refluxing at patch boundaries
  // FineAdd: If we are not on the coarsest level, the fluxes at this level will form part of the correction
  //          to a coarse level
  // CrseInit:  If we are not the finest level, the fluxes at patch boundaries need correcting.  Since we
  //            know that the coarse level happens first, we initialise the boundary fluxes through this
  //            function, and subsequently FineAdd will modify things ready for the correction
  // Both functions have the same arguments:
  // First: Name of the flux MultiFab (this is done dimension-by-dimension)
  // Second: Direction, to ensure the correct vertices are being corrected
  // Third: Source component - the first entry of the flux MultiFab that is to be copied (it is possible that
  //        some variables will not need refluxing, or will be computed elsewhere
  // Fourth: Destination component - the first entry of the flux register that this call to FineAdd sends to
  // Fifth: NUM_STATE - number of states being added to the flux register
  // Sixth: Multiplier - in general, the least accurate (coarsest) flux is subtracted (-1) and the most accurate (finest) flux is added (+1)
  if (do_reflux)
  {
    if (current)
    {
      for (int i = 0; i < amrex::SpaceDim; i++)
        current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, 1.);
    }
    if (fine)
    {
      for (int i = 0; i < amrex::SpaceDim; i++)
        fine->CrseInit(fluxes[i], i, 0, 0, NUM_STATE, -1.);
    }
  }

  return dt;
}

// Estimate time step.
// This function is called by all of the other time step functions in AMReX
Real AmrLevelAdv::estTimeStep(Real)
{
  // This is just a dummy value to start with
  Real dt_est = 1.0e+20;
  State uVec, wVec;

  const Real *dx = geom.CellSize();
  const Real *prob_lo = geom.ProbLo();
  const Real cur_time = state[Phi_Type].curTime();
  MultiFab &S_new = get_new_data(Phi_Type);

  Real maxSpeed, vx, vy, vz, v;
  
  for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {
      const Box &bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto &arr = S_new.array(mfi);

      for (int k = lo.z; k <= hi.z; k++)
      {
        for (int j = lo.y; j <= hi.y; j++)
        {
          for (int i = lo.x; i <= hi.x; i++)
          {
            for (int ii = 0; ii < NUM_STATE; ii++) {
              uVec[ii] = arr(i, j, k, ii);
            }
            maxSpeed = std::max(maxSpeed, fastSpeed(uVec, 0) + std::abs(uVec[MOM_X]) / uVec[RHO]);
            maxSpeed = std::max(maxSpeed, fastSpeed(uVec, 1) + std::abs(uVec[MOM_Y]) / uVec[RHO]);
            maxSpeed = std::max(maxSpeed, fastSpeed(uVec, 2) + std::abs(uVec[MOM_Z]) / uVec[RHO]);
          }
        }
      }
    }

  for (unsigned int d = 0; d < amrex::SpaceDim; ++d)
  {
    dt_est = std::min(dt_est, dx[d] / maxSpeed);
  }

  // Ensure that we really do have the minimum across all processors
  ParallelDescriptor::ReduceRealMin(dt_est);
  dt_est *= cfl;

  if (verbose)
  {
    amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level << ":  dt_est = " << dt_est << std::endl;
  }

  return dt_est;
}

// Compute initial time step.
Real AmrLevelAdv::initialTimeStep()
{
  return estTimeStep(0.0);
}

// Compute initial dt
void AmrLevelAdv::computeInitialDt(int finest_level, int sub_cycle, Vector<int> &n_cycle, const Vector<IntVect> &ref_ratio, Vector<Real> &dt_level, Real stop_time)
{
  // Grids have been constructed, compute dt for all levels.
  if (level > 0)
    return;

  // Initial guess
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0, n_factor * dt_level[i]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 0.001 * dt_0;
  Real cur_time = state[Phi_Type].curTime();
  if (stop_time >= 0.0)
  {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

// Compute new dt
void AmrLevelAdv::computeNewDt(int finest_level, int sub_cycle, Vector<int> &n_cycle, const Vector<IntVect> &ref_ratio, Vector<Real> &dt_min, Vector<Real> &dt_level, Real stop_time, int post_regrid_flag)
{
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  if (level > 0)
    return;

  // Although we only compute the time step on the finest level, we
  // need to take information from all levels into account.  The
  // sharpest features may be smeared out on coarse levels, so not
  // using finer levels could cause instability
  for (int i = 0; i <= finest_level; i++)
  {
    AmrLevelAdv &adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  // A couple of things are implemented to ensure that timesteps
  // don't suddenly grow by a lot, as this could lead to errors
  if (post_regrid_flag == 1)
  {
    // Limit dt's by pre-regrid dt
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i], dt_level[i]);
    }
  }
  else
  {
    // Limit dt's by change_max * old dt
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
    }
  }

  // Find the minimum over all levels
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
  }

  // Limit dt's by the value of stop_time
  const Real eps = 0.001 * dt_0;
  Real cur_time = state[Phi_Type].curTime();
  if (stop_time >= 0.0)
  {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

// Do work after timestep().
// If something has to wait until all processors have done their advance function, the post_timestep function
// is the place to put it.  Refluxing and averaging down are two standard examples for AMR
void AmrLevelAdv::post_timestep(int iteration)
{
  // Integration cycle on fine level grids is complete
  // do post_timestep stuff here.
  int finest_level = parent->finestLevel();

  if (do_reflux && level < finest_level)
    reflux();

  if (level < finest_level)
    avgDown();
}

// Do work after regrid().
void AmrLevelAdv::post_regrid(int lbase, int new_finest)
{
}

// Do work after a restart().
void AmrLevelAdv::post_restart()
{
}

// Do work after init().
void AmrLevelAdv::post_init(Real stop_time)
{
  if (level > 0)
    return;
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  int finest_level = parent->finestLevel();
  for (int k = finest_level - 1; k >= 0; k--)
    getLevel(k).avgDown();
}

// Error estimation for regridding.
//  Determine which parts of the domain need refinement
void AmrLevelAdv::errorEst(TagBoxArray &tags, int clearval, int tagval, Real time, int n_error_buf, int ngrow)
{
  const Real *dx = geom.CellSize();
  const Real *prob_lo = geom.ProbLo();

  MultiFab &S_new = get_new_data(Phi_Type);
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  Sborder.FillBoundary(geom.periodicity());

  Vector<BCRec> bc(NUM_STATE); // copied from AmReX documentation
  for (int n = 0; n < NUM_STATE; ++n)
  {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      {
          if (geom.isPeriodic(idim))
          {
              bc[n].setLo(idim, BCType::int_dir); // interior
              bc[n].setHi(idim, BCType::int_dir);
          }
          else
          {
              bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation
              bc[n].setHi(idim, BCType::foextrap);
          }
      }
  }

  FillDomainBoundary(Sborder, geom, bc);
  Vector<int> itags;

  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
  {
    const Box &tilebx = mfi.tilebox();

    // An AMReX construction, effectively a boolean array which is true in positions that are valid for refinement
    TagBox &tagfab = tags[mfi];

    tagfab.get_itags(itags, tilebx);

    // data pointer and index space
    int *tptr = itags.dataPtr();
    const int *tlo = tilebx.loVect();
    const int *thi = tilebx.hiVect();

    state_error(tptr, AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
                BL_TO_FORTRAN_3D(Sborder[mfi]),
                &tagval, &clearval,
                AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()),
                AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);

    // Now update the tags in the TagBox.
    tagfab.tags_and_untags(itags, tilebx);
  }
}

// This function reads the settings file
void AmrLevelAdv::read_params()
{
  // Make sure that this is only done once
  static bool done = false;

  if (done)
    return;

  done = true;

  // A ParmParse object allows settings, with the correct prefix, to be read in from the settings file
  // The prefix can help identify what a settings parameter is used for
  // AMReX has some default ParmParse names, amr and geometry are two commonly needed ones
  ParmParse pp("adv");

  // ParmParse has two options; query and get.  Query will only alter
  // a parameter if it can be found (if these aren't in the settings
  // file, then the values at the top of this file will be used).  Get
  // will throw an error if the parameter is not found in the settings
  // file.
  pp.query("v", verbose);
  pp.query("cfl", cfl);
  pp.query("do_reflux", do_reflux);

  Vector<Real> wL0Help, wR0Help;
  pp.get("Gamma", Gamma);
  pp.getarr("wL0", wL0Help);
  pp.getarr("wR0", wR0Help);
  for (int i = 0; i < NUM_STATE; i++) {
    wL0[i] = wL0Help[i];
    wR0[i] = wR0Help[i];
  }
  pp.get("orientation", orient);
  if (order == 1) {
    limFunc = forcezero; // no reconstruction performed in the first order case
  }
  else {
    limFunc = minbee;
  }
  std::string limiting;
  pp.get("limiting", limiting);
  if (limiting == "E") {
    for (int i = 0; i < NUM_STATE; i++) {
      varLim[i] = ENE;
    }
  }
  else {
    for (int i = 0; i < NUM_STATE; i++) {
      varLim[i] = i;
    }
  }

  Geometry const *gg = AMReX::top()->getDefaultGeometry();

  if (!gg->IsCartesian())
  {
    amrex::Abort("Please set geom.coord_sys = 0");
  }

  // read tagging parameters from probin file
  // Tradtionally, the inputs file with ParmParse functionality is handled by C++.  However, a Fortran settings
  // file, by default named probin, can also supply variables.  Mostly used for mesh refinement (tagging) criteria
  std::string probin_file("probin");

  ParmParse ppa("amr");
  ppa.query("probin_file", probin_file);

  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  // use a fortran routine to
  // read in tagging parameters from probin file
  get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);
}

// AMReX has an inbuilt reflux command, but we still have the freedom
// to decide what goes into it (for example, which variables are
// actually refluxed).  This also gives a little flexibility as to
// where flux registers are stored
void AmrLevelAdv::reflux()
{
  BL_ASSERT(level < parent->finestLevel());

  const Real strt = amrex::second();

  // Call the reflux command with the appropriate data.  Because there
  // are no flux registers on the coarse level, they start from the
  // first level.  But the coarse level to the (n-1)^th are the ones
  // that need refluxing, hence the `level+1'.
  getFluxReg(level + 1).Reflux(get_new_data(Phi_Type), 1.0, 0, 0, NUM_STATE, geom);

  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real end = amrex::second() - strt;

    ParallelDescriptor::ReduceRealMax(end, IOProc);

    amrex::Print() << "AmrLevelAdv::reflux() at level " << level << " : time = " << end << std::endl;
  }
}

// Generic function for averaging down - in this case it just makes sure it doesn't happen on the finest level
void AmrLevelAdv::avgDown()
{
  if (level == parent->finestLevel())
  {
    return;
  }
  // Can select which variables averaging down will happen on - only one to choose from in this case!
  avgDown(Phi_Type);
}

// Setting up the call to the AMReX-implemented average down function
void AmrLevelAdv::avgDown(int state_indx)
{
  // For safety, again make sure this only happens if a finer level exists
  if (level == parent->finestLevel())
    return;

  AmrLevelAdv &fine_lev = getLevel(level + 1);
  MultiFab &S_fine = fine_lev.get_new_data(state_indx);
  MultiFab &S_crse = get_new_data(state_indx);

  // Call the AMReX average down function:
  // S_fine: Multifab with the fine data to be averaged down
  // S_crse: Multifab with the coarse data to receive the fine data where necessary
  // fine_lev.geom:  Geometric information (cell size etc.) for the fine level
  // geom: Geometric information for the coarse level (i.e. this level)
  // 0: First variable to be averaged (as not all variables need averaging down
  // S_fine.nComp(): Number of variables to average - this can be computed automatically from a multifab
  // refRatio: The refinement ratio between this level and the finer level
  amrex::average_down(S_fine, S_crse, fine_lev.geom, geom, 0, S_fine.nComp(), parent->refRatio(level));
}
