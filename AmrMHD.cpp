#include "AmrMHD.H"

using namespace amrex;

typedef FixedArray<Real, 3> Vec3;

// Functions for solving ideal MHD equations
Real BSquared(State uwVec) {
  // Returns the square of the magnetic field magnitude. uwVec can be the vector of primitive or conservative variables
  return uwVec[BX] * uwVec[BX] + uwVec[BY] * uwVec[BY] + uwVec[BZ] * uwVec[BZ];
}

Real energy(State wVec) {
  // Calculates total energy based on vector of primitive variables
  Real vx = wVec[VX], vy = wVec[VY], vz = wVec[VZ];
  Real v2 = vx * vx + vy * vy + vz * vz;
  Real B2 = BSquared(wVec);
  return wVec[PRE] / (Gamma - 1) + 0.5 * wVec[RHO] * (v2 + B2);
}

Real pressure(State uVec) {
  // Calculates pressure based on vector of conservative variables
  Real rho = uVec[RHO], mom_x = uVec[MOM_X], mom_y = uVec[MOM_Y], mom_z = uVec[MOM_Z];
  Real B2 = BSquared(uVec);
  Real v2 = (mom_x * mom_x + mom_y * mom_y + mom_z * mom_z) / rho;
  return (Gamma - 1) * (uVec[ENE] - 0.5 * (v2 + B2));
}

void flux(State uVec, State& fluxVec, unsigned int coord) { // gotta change this for 3D
  // Calculates flux along a given axis
  Real p = pressure(uVec);
  Real rho = uVec[RHO], mom_x = uVec[MOM_X], mom_y = uVec[MOM_Y], mom_z = uVec[MOM_Z], E = uVec[ENE], Bx = uVec[BX], By = uVec[BY], Bz = uVec[BZ];
  Real B2 = BSquared(uVec);
  Real pT = p + 0.5 * B2;
  unsigned int vCoord = 1 + coord, BCoord = 5 + coord;
  Real mom_Norm = uVec[vCoord];
  Real BNorm = uVec[BCoord];
  fluxVec[0] = mom_Norm;
  fluxVec[1] = mom_x * mom_Norm / rho - Bx * BNorm;
  fluxVec[2] = mom_y * mom_Norm / rho - By * BNorm;
  fluxVec[vCoord] += pT;
  fluxVec[3] = E * mom_Norm / rho - Bz * BNorm;
  fluxVec[4] = (mom_Norm * (E + pT) - BNorm * (mom_x * Bx + mom_y * By + mom_z * Bz)) / rho;
  fluxVec[BCoord] = 0.;
  fluxVec[6 - coord] = (mom_x * By - mom_y * Bx) / rho * pow(-1, coord);
  fluxVec[7] = (mom_Norm * Bz - BNorm * u4) / rho;
  // fluxVec[8] = 0.;
}

void primitive(State uVec, State& wVec) {
  // Converts from conservative to primitive variables
  Real rho = uVec[RHO];
  wVec = uVec;
  wVec[VX] /= rho;
  wVec[VY] /= rho;
  wVec[VZ] /= rho;
  wVec[PRE] = pressure(uVec);
}

void conservative(State wVec, State& uVec) {
  // Converts from primitive to conservative variables
  Real rho = wVec[RHO];
  uVec = wVec;
  uVec[MOM_X] *= rho;
  uVec[MOM_Y] *= rho;
  uVec[MOM_Z] *= rho;
  uVec[ENE] = energy(wVec);
}

Real fastSpeed(State uVec, unsigned int coord) {
  // Calculates the magnetoacoustic fast speed based on conservative vector of variables
  Real p = pressure(uVec), rho = uVec[RHO], Bx = uVec[BX], By = uVec[BY], Bz = uVec[BZ];
  Real B2 = BSquared(uVec);
  Real c_a2 = B2 / rho;
  Real c_s2 = Gamma * p / rho;
  Real speedSum = c_a2 + c_s2;
  Real BNorm = uVec[5 + coord];
  return sqrt(0.5 * (speedSum + sqrt(speedSum * speedSum - 4 * c_s2 * BNorm * BNorm / rho)));
}

void estimateWaveSpeeds(State wL, State wR, Vec3& waveSpeeds, unsigned int coord) {
  /* Estimates HLLC wave speeds */
  Real SL, SR, SStar, maxFast;
  unsigned int vCoord = 1 + coord, BCoord = 5 + coord;
  Real rhoL = wL[RHO], vNormL = wL[vCoord], pL = wL[PRE], BxL = wL[BX], ByL = wL[BY], BzL = wL[BZ]; 
  Real rhoR = wR[RHO], vNormR = wR[vCoord], pR = wR[PRE], BxR = wR[BX], ByR = wR[BY], BzR = wR[BZ];
  Real BNormR = wR[BCoord], BNormL = wL[BCoord];
  Real B2R = BSquared(wR);
  Real B2L = BSquared(wL);
  Real pModL = pL + 0.5 * B2L - BNormL * BNormL;
  Real pModR = pR + 0.5 * B2R - BNormR * BNormR;
  maxFast = std::max(fastSpeed(wL, coord), fastSpeed(wR, coord));
  SL = std::min(vNormL, vNormR) - maxFast;
  SR = std::max(vNormL, vNormR) + maxFast;
  SStar = (rhoR * vNormR * (SR - vNormR) - rhoL * vNormL * (SL - vNormL) + pModL - pModR) / (rhoR * (SR - vNormR) - rhoL * (SL - vNormL));
  waveSpeeds[0] = SL;
  waveSpeeds[1] = SStar;
  waveSpeeds[2] = SR;
}

void calcUStarK(State wK, State wHLL, Real EK, Real SK, Real SStar, State& uStarK, unsigned int coord) {
  /* Calculates HLLC intermediate state */
  Real rhoK = wK[RHO], vxK = wK[VX], vyK = wK[VY], vzK = wK[VZ], pK = wK[PRE], BxK = wK[BX], ByK = wK[BY], BzK = wK[BZ];
  Real vNormK = wK[1 + coord], BNormK = wK[5 + coord];
  Real rhoStarK = rhoK * (SK - vNormK) / (SK - SStar);
  Real BNormStar = wK[5 + coord];
  Real pKT = pK + 0.5 * BSquared(wK);
  Real pStarT = rhoK * (SK - vNormK) * (SStar - vNormK) + BNormStar * BNormStar + pKT - BNormK * BNormK;
  Real vdotBK = vxK * BxK + vyK * ByK + vzK * BzK;
  Real vdotBHLL = wHLL[VX] * wHLL[BX] + wHLL[VY] * wHLL[BY] + wHLL[VZ] * wHLL[BZ];
  uStarK[0] = rhoStarK;
  uStarK[1 + coord] = rhoStarK * SStar;
  uStarK[2 - coord] = rhoStarK * wK[2 - coord] - (BNormStar * wHLL[6 - coord] - BxK * ByK) / (SK - SStar);
  uStarK[3] = rhoStarK * vzK - (BNormStar * wHLL[BZ] - BNormK * BzK) / (SK - SStar);
  uStarK[4] = EK * (rhoStarK / rhoK) + ((pStarT * SStar - pKT * vNormK) - (BNormStar * vdotBHLL - BNormK * vdotBK)) / (SK - SStar);
  uStarK[5 + coord] = BNormStar;
  uStarK[6 - coord] = wHLL[6 - coord];
  uStarK[7] = wHLL[BZ];
}

void HLLCFlux(State uL, State uR, State& fluxVec, unsigned int coord) {
  /* Calculates HLLC flux along a given axis */
  Real SL, SStar, SR;
  State wL, wR, uStarK, fluxL, fluxR, uHLL, wHLL;
  Vec3 waveSpeeds;
  primitive(uL, wL);
  primitive(uR, wR);
  estimateWaveSpeeds(wL, wR, waveSpeeds, coord);
  SL = waveSpeeds[0], SStar = waveSpeeds[1], SR = waveSpeeds[2];
  if (SL > 0.) {
    flux(uL, fluxVec, coord);
  }
  else if ((SStar >= 0) || (SStar <= 0. && SR >= 0.)) {
    flux(uL, fluxL, coord);
    flux(uR, fluxR, coord);
    for (int i = 0; i < NUM_STATE; i++) {
        uHLL[i] = (SR * uR[i] - SL * uL[i] - fluxR[i] + fluxL[i]) / (SR - SL);
    }
    primitive(uHLL, wHLL);
    State wK, uK, fluxK;
    Real SK;
    if (SStar >= 0) {
        wK = wL, uK = uL, SK = SL, fluxK = fluxL;
    }
    else {
        wK = wR, uK = uR, SK = SR, fluxK = fluxR;
    }
    calcUStarK(wK, wHLL, uK[ENE], SK, SStar, uStarK, coord);
    for (int i = 0; i < NUM_STATE; i++) {
        fluxVec[i] = fluxK[i] + SK * (uStarK[i] - uK[i]);
    }
  }
  else {
    flux(uR, fluxVec, coord);
  }
}

void reconstruct(State uL, State u0, State uR, State& u0Re_L, State& u0Re_R, State varLim) {
  /* Reconstruct state u0 at boundaries. uL/uR = state at the cell to the left/right of u0, u0Re_L/R = reconstructed state u0 at the left/right interfaces */
  Real r, uDelta_i, eps = 1e-9, limDeltaL, limDeltaR, xi;
  for (int i = 0; i < NUM_STATE; i++) {
    limDeltaR = uR[varLim[i]] - u0[varLim[i]];
    if (std::abs(limDeltaR) < eps) {
        u0Re_L[i] = u0Re_R[i] = u0[i];
    }
    else {
        limDeltaL = u0[varLim[i]] - uL[varLim[i]];
        r = limDeltaL / limDeltaR;
        xi = limFunc(r);
        uDelta_i = 0.25 * xi * (uR[i] - uL[i]);
        u0Re_L[i] = u0[i] - uDelta_i;
        u0Re_R[i] = u0[i] + uDelta_i;
    }
  }
}

void halfTimeEvol(State& u0Re_L, State& u0Re_R, Real step, unsigned int coord) {
  /* Performs half time local evolution according to the MUSCL-Hancock scheme */
  State fluxR, fluxL;
  Real correction;
  flux(u0Re_L, fluxL, coord);
  flux(u0Re_R, fluxR, coord);
  for (int i = 0; i < NUM_STATE; i++) {
    correction = 0.5 * step * (fluxR[i] - fluxL[i]);
    u0Re_L[i] -= correction;
    u0Re_R[i] -= correction;
  }
}