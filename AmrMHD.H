#ifndef AmrMHD_H
#define AmrMHD_H
#include "AMReX_REAL.H"
#include "AMReX_Array.H"

using namespace amrex;

enum ConsVar{RHO=0, MOM_X, MOM_Y, MOM_Z, ENE, BX, BY, BZ};
enum PrimVar{VX=1, VY, VZ, PRE};

static const int stateVars = 8;

typedef GpuArray<int, stateVars> varMap;
typedef GpuArray<Real, stateVars> State;
typedef GpuArray<Real, 3> Vec3;

Real BSquared(State uwVec);
Real energy(State wVec);
Real pressure(State uVec);
void flux(State uVec, State& fluxVec, unsigned int coord);
void primitive(State uVec, State& wVec);
void conservative(State wVec, State& uVec);
Real fastSpeed(State uVec, unsigned int coord);
void estimateWaveSpeeds(State wL, State wR, Vec3& waveSpeeds, unsigned int coord);
void calcUStarK(State wK, State wHLL, Real EK, Real SK, Real SStar, State& uStarK, unsigned int coord);
void HLLCFlux(State uL, State uR, State& fluxVec, unsigned int coord);
void reconstruct(State uL, State u0, State uR, State& u0Re_L, State& u0Re_R, varMap varLim);
void halfTimeEvol(State& u0Re_L, State& u0Re_R, Real step, unsigned int coord);

#endif