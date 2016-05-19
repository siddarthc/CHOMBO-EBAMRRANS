#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchTransport.H"
#include "EBPatchTransportF_F.H"
#include "EBLoHiCenter.H"
#include "EBArith.H"
#include "PolyGeom.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "CH_Timer.H"
#include "EBPatchGodunovF_F.H"
#include "EBScalarAdvectBC.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

/******/
EBPatchTransportFactory::
EBPatchTransportFactory(RefCountedPtr<EBPhysIBCFactory>& a_bcFactoryPtr,
                     const bool&                      a_useLimiting,
                     const Real&                      a_maxVal,
                     const Real&                      a_minVal,
                     const bool&                      a_setMaxMin,
                     const bool&                      a_useFourthOrderSlopes,
                     const bool&                      a_useFlattening,
                     const bool&                      a_useArtificialVisc)
  :EBPatchGodunovFactory()
{
  m_bcFactoryPtr = a_bcFactoryPtr;
  m_useLimiting  = a_useLimiting;
  m_setMaxMin    = a_setMaxMin;
  m_maxVal       = a_maxVal;
  m_minVal       = a_minVal;
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_useFlattening = a_useFlattening;
  m_useArtificialVisc = a_useArtificialVisc;
}
/******************/
EBPatchTransportFactory::
~EBPatchTransportFactory()
{
}
/******************/
EBPatchGodunov*
EBPatchTransportFactory::
create() const
{
  EBPatchTransport* retval= new EBPatchTransport();
  retval->setEBPhysIBC(*m_bcFactoryPtr);
  retval->useLimiting(m_useLimiting);
  retval->setSlopeParameters(m_useFourthOrderSlopes,
                             m_useFlattening,
                             m_useLimiting);
  retval->artificialViscosity(m_useArtificialVisc);
  if (m_setMaxMin)
    {
      retval->setMaxMin(m_maxVal, m_minVal);
    }
  return static_cast<EBPatchGodunov*>(retval);
}
/******/
void
EBPatchTransport::
setValidBox(const Box&        a_validBox,
            const EBISBox&    a_ebisBox,
            const IntVectSet& a_coarseFineIVS,
            const Real&       a_cumulativeTime,
            const Real&       a_timeStep)
{
  CH_TIME("EBPatchTransport::setValidBox");
  EBPatchGodunov::setValidBox(a_validBox, a_ebisBox, a_coarseFineIVS, a_cumulativeTime, a_timeStep);
  Box slopeBoxG2 = grow(a_validBox, 2);
  slopeBoxG2 &= m_domain;
  IntVectSet ivs = m_ebisBox.getIrregIVS(slopeBoxG2);
  m_vofsStencil.define(ivs, m_ebisBox.getEBGraph(), 1);
  for(VoFIterator vofit(ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      int radius = 1;
      EBArith::getAllVoFsInMonotonePath(m_vofsStencil(vofit(), 0), vofit(), m_ebisBox, radius);
    }
}
/******************/
bool EBPatchTransport::
checkCoveredFaceStencil(const VolIndex& a_vof,
                        const Vector<VolIndex>& a_vofsStencil,
                        const RealVect& a_normal,
                        const int&      a_faceDir,
                        const Side::LoHiSide& a_sd)
{
  bool hasAllVoFs = false;
  if(SpaceDim==2)
    {
      hasAllVoFs = checkCoveredFaceStencil2d(a_vof, a_vofsStencil, a_normal, a_faceDir, a_sd);
    }
  else
    {
      hasAllVoFs = checkCoveredFaceStencil3d(a_vof, a_vofsStencil, a_normal, a_faceDir, a_sd);
    }
  return hasAllVoFs;
}
/******************/
bool EBPatchTransport::
checkCoveredFaceStencil2d(const VolIndex& a_vof,
                          const Vector<VolIndex>& a_vofsStencil,
                          const RealVect& a_normal,
                          const int&      a_faceDir,
                          const Side::LoHiSide& a_sd)
{
  CH_assert(SpaceDim==2);
  bool hasAllVoFs;
  int tangenDir = 1 - a_faceDir;

  int signNorm = 1;
  int signTang = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  if (Abs(a_normal[tangenDir]) == 0.0)
    {//since signTang is arbitrary for this case,
      //check if we are next to the high domain edge
      IntVect ivHi = a_vof.gridIndex();;
      ivHi[tangenDir] += 1;
      if (!m_domain.contains(ivHi))
        {
          signTang = -1;
        }
    }
  else if (a_normal[tangenDir] < 0.0)
    {
      signTang = -1;
    }

  //iv[0][0] is the vof. iv[1][0] is vofside.  iv[0][1] is vofup iv[1][1] is vofcorner
  IntVect   ivSten[2][2];
  bool      hasVoF[2][2];
  VolIndex     vof[2][2];
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], a_vofsStencil, ivSten[ix][iy]);
        }
    }

  int dirBigNorm, dirLitNorm;
  bool hasVoFBigNorm,  hasVoFLitNorm;
  int  signBigNorm, signLitNorm;
  Real eps = 1.0e-12;
  if ((Abs(a_normal[a_faceDir]) < eps) || (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir]))
    {
      dirBigNorm = tangenDir;
      dirLitNorm = a_faceDir;

      signBigNorm   = signTang;
      signLitNorm   = signNorm;
      hasVoFBigNorm = hasVoF[0][1];
    }
  else
    {
      dirBigNorm = a_faceDir;
      dirLitNorm = tangenDir;

      signBigNorm   = signNorm;
      signLitNorm   = signTang;
      hasVoFBigNorm = hasVoF[1][0];

    }
  hasVoFLitNorm = hasVoF[1][1];
  hasAllVoFs = (hasVoFLitNorm && hasVoFBigNorm);

  return hasAllVoFs;
}
/*********/
bool EBPatchTransport::
checkCoveredFaceStencil3d(const VolIndex& a_vof,
                          const Vector<VolIndex>& a_vofsStencil,
                          const RealVect& a_normal,
                          const int&      a_faceDir,
                          const Side::LoHiSide& a_sd)
{
  CH_assert(SpaceDim==3);
  bool hasAllVoFs;
  Tuple<int,CH_SPACEDIM-1> tangenDir = PolyGeom::computeTanDirs(a_faceDir);
  Real anormNorm = Abs(a_normal[a_faceDir]);
  Real anormTan[2];
  int signNorm = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  int signTang[2];
  for (int itan = 0; itan < 2; itan++)
    {
      int tanDir = tangenDir[itan];
      anormTan[itan] = Abs(a_normal[tanDir]);
      signTang[itan] =  1;
      if (anormTan[itan] == 0.0)
        {//since signTang is arbitrary for this case,
          //check if we are next to the high domain edge
          IntVect ivHi = a_vof.gridIndex();;
          ivHi[tanDir] += 1;
          if (!m_domain.contains(ivHi))
            {
              signTang[itan] = -1;
            }
        }
      else if (a_normal[tanDir] < 0.0)
        {
          signTang[itan] = -1;
        }
    }

  const IntVect& ivVoF= a_vof.gridIndex();

  //whice one of the tangential directions has the largest normal
  int d1, d2;
  if (anormTan[0]/m_dx[tangenDir[0]] > anormTan[1]/m_dx[tangenDir[1]])
    {
      d1 = 0;
      d2 = 1;
    }
  else
    {
      d1 = 1;
      d2 = 0;
    }

  // figure out in which plane we are extrapolating
  bool faceDirOut = ((anormNorm/m_dx[a_faceDir] > anormTan[0]/m_dx[tangenDir[0]]) &&
                     (anormNorm/m_dx[a_faceDir] > anormTan[1]/m_dx[tangenDir[1]]));

  IntVect ivSten[2][2];
  bool hasVoF[2][2];
  VolIndex vofSten[2][2];

  if (faceDirOut)
    {// face direction has largest normal.

      ivSten[0][0] = ivVoF + signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[0]*BASISV(tangenDir[0]);
      ivSten[0][1] = ivVoF + signTang[1]*BASISV(tangenDir[1]);
      ivSten[1][1] = ivVoF + signTang[0]*BASISV(tangenDir[0]) + signTang[1]*BASISV(tangenDir[1]) ;
      d1 = 0;
      d2 = 1;

    }
  else
    { //tandir[d1] is the biggest

      ivSten[0][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      - signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      ;
      ivSten[0][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) - signNorm*BASISV(a_faceDir);
      ivSten[1][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) ;

    }

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], a_vofsStencil, ivSten[ix][iy]);
        }
    }
  hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  return hasAllVoFs;
}
/******************/
Real
EBPatchTransport::
artificialViscosityCoefficient() const
{
/*
  ParmParse pp;
  Real retval;
  pp.get("scal_artificial_viscosity", retval);
  return retval;
*/
  MayDay::Error("should not be called");
  return -1.0;
}
/******/
void EBPatchTransport::
normalPred(EBCellFAB&       a_rhoLo,
           EBCellFAB&       a_rhoHi,
           const EBCellFAB& a_rho,
           const EBCellFAB& a_dRho,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_assert(m_isNormalVelSet);
  CH_assert(m_normalVelPtr != NULL);
  const EBCellFAB& velcc = *m_normalVelPtr;

  int ivar = 0;

  /**/
  BaseFab<Real>& regRhoLo    = a_rhoLo.getSingleValuedFAB();
  BaseFab<Real>& regRhoHi    = a_rhoHi.getSingleValuedFAB();
  const BaseFab<Real>& regVel =  velcc.getSingleValuedFAB();
  const BaseFab<Real>& regRho =  a_rho.getSingleValuedFAB();
  const BaseFab<Real>& regDRho= a_dRho.getSingleValuedFAB();

  FORT_PREDADVECT(CHF_BOX(a_box),
                  CHF_CONST_FRA1(regRho,ivar),
                  CHF_CONST_FRA1(regDRho, ivar),
                  CHF_CONST_FRA(regVel),
                  CHF_FRA1(regRhoLo, ivar),
                  CHF_FRA1(regRhoHi, ivar),
                  CHF_CONST_INT(a_dir),
                  CHF_CONST_REAL(a_dtbydx));
  /**/

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real dense, denlo, denhi, denslope;
      RealVect veloc;

      dense    = a_rho(vof, ivar);
      denslope = a_dRho(vof, ivar);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = velcc(vof, idir);
        }

      FORT_POINTPREDADVECT(CHF_REAL(dense),
                           CHF_REAL(denlo),
                           CHF_REAL(denhi),
                           CHF_REAL(denslope),
                           CHF_REALVECT(veloc),
                           CHF_CONST_INT(a_dir),
                           CHF_CONST_REAL(a_dtbydx));

      a_rhoLo(vof, ivar) = denlo;
      a_rhoHi(vof, ivar) = denhi;
    } 
}
/*****************************/
bool
EBPatchTransport::
usesFlattening() const
{
//  CH_assert(m_isSlopeSet);
//  return m_useFlattening;
//  return always false for now
  return false;
}
/******/
bool
EBPatchTransport::
usesArtificialViscosity() const
{
  return false;
}
/******/
bool
EBPatchTransport::
usesFourthOrderSlopes() const
{
  CH_assert(m_isSlopeSet);
  return m_useFourthOrderSlopes;
}
/******/
void
EBPatchTransport::
consToPrim(EBCellFAB&       a_primState,
           const EBCellFAB& a_consState,
           const Box&       a_box,
           int              a_logflag,
           bool             a_verbose)
{
  CH_TIME("EBPatchTransport::consToPrim");
  CH_assert(numPrimitives() == numConserved());
  CH_assert(a_primState.nComp() == numPrimitives());
  CH_assert(a_primState.nComp() == a_consState.nComp());
  CH_assert(a_primState.getRegion().contains(a_box));
  Interval interv(0, numConserved()-1);
  a_primState.copy(a_box, interv, a_box, a_consState, interv);
}
/******/
void
EBPatchTransport::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const EBCellFAB&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchTransport::consToPrimIrr");
  CH_assert(numPrimitives() == numConserved());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_primState(vofit(), ivar) = a_consState(vofit(), ivar);
        }
    }
}
/******/
void
EBPatchTransport::
consToPrim(BaseIVFAB<Real>&        a_primState,
           const BaseIVFAB<Real>&  a_consState,
           const IntVectSet&       a_ivs)
{
  CH_TIME("EBPatchTransport::consToPrimIrrIrr");
  CH_assert(numPrimitives() == numConserved());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_primState(vofit(), ivar) = a_consState(vofit(), ivar);
        }
    }
}
/******/
void
EBPatchTransport::
primToCons(EBCellFAB&       a_consState,
           const EBCellFAB& a_primState,
           const Box&       a_box)
{
  CH_TIME("EBPatchTransport::primToConsRegReg");
  CH_assert(numPrimitives() == numConserved());
  Interval interv(0, numConserved()-1);
  a_consState.copy(a_box, interv, a_box, a_primState, interv);
}
/******/
void
EBPatchTransport::
primToCons(BaseIVFAB<Real>&       a_consState,
           const BaseIVFAB<Real>& a_primState,
           const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchTransport::primToConsRegIrr");
  CH_assert(numPrimitives() == numConserved());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_consState(vofit(), ivar) = a_primState(vofit(), ivar);
        }
    }
}
/***/
void
EBPatchTransport::
riemann(EBFaceFAB&       a_primGdnv,
        const EBCellFAB& a_primLeft,
        const EBCellFAB& a_primRigh,
        const int&       a_faceDir,
        const Box&       a_box)
{
  //these two need to be set
  CH_assert((s_doingVel == 1) || (s_doingVel == 0));
  //  CH_assert((s_curComp >= 0) && (s_curComp < SpaceDim);
  CH_assert(s_curComp >= 0);
  CH_assert(m_isAdvVelSet);
  CH_assert(!m_domain.isPeriodic(a_faceDir));

  // const EBCellFAB& normalVel = *m_normalVelPtr;
  const EBFluxFAB& advectVel = *m_advVelPtr;
  const EBFaceFAB& advectVelFaceFAB = advectVel[a_faceDir];
  int nPrim = numPrimitives();

  Box cellBox = enclosedCells(a_box);

  //explicitly cast the left and right states to modifiable references
  //because we need to shift them to faces and then shift them back.
  BaseFab<Real>& regPrimRigh = (BaseFab<Real>&)a_primRigh.getSingleValuedFAB();
  BaseFab<Real>& regPrimLeft = (BaseFab<Real>&)a_primLeft.getSingleValuedFAB();
  BaseFab<Real>& regPrimGdnv = a_primGdnv.getSingleValuedFAB();
  // const BaseFab<Real>& regNormalVel = normalVel.getSingleValuedFAB();
  const BaseFab<Real>& regAdvectVel = advectVelFaceFAB.getSingleValuedFAB();
  //find the regular part of qgdnv.
  FORT_ADVECTRIEMANN(CHF_BOX(a_box),
                     CHF_FRA(regPrimGdnv),
                     CHF_CONST_FRA(regPrimLeft),
                     CHF_CONST_FRA(regPrimRigh),
                     // CHF_CONST_FRA(regNormalVel),
                     CHF_CONST_FRA1(regAdvectVel, 0),//only 1 (normal) component at face
                     CHF_CONST_INT(a_faceDir),
                     CHF_CONST_INT(nPrim),
                     CHF_CONST_INT(s_curComp),
                     CHF_CONST_INT(s_doingVel));

  //the box sent into this is face-centered.
  //we need to use the cell-centered one it surrounds.
  Vector<FaceIndex> multiFaces =  m_ebisBox.getEBGraph().getMultiValuedFaces(a_faceDir, cellBox);
  for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (!face.isBoundary())
        {
          VolIndex vofl = face.getVoF(Side::Lo);
          VolIndex vofr = face.getVoF(Side::Hi);

          Real velFace = advectVelFaceFAB(face, 0);
          for (int ivar = 0; ivar < nPrim; ivar++)
            {
              if (velFace > TOL)
                {
                  a_primGdnv(face, ivar) = a_primLeft(vofl, ivar);
                }
              else if (velFace < -TOL)
                {
                  a_primGdnv(face, ivar) = a_primRigh(vofr, ivar);
                }
              else
                {
                  if ((s_doingVel == 1) && (s_curComp == a_faceDir))
                    {
                      a_primGdnv(face, ivar) = 0.0;
                    }
                  else
                    {
                      a_primGdnv(face, ivar) = 0.5*(a_primRigh(vofr, ivar) + a_primLeft(vofl, ivar));
                    }
                }
            }

        }
    }
}
/***/
/*****************************/
void
EBPatchTransport::
riemann(BaseIVFAB<Real>&        a_coveredPrim,
        const BaseIVFAB<Real>&  a_exteState,
        const EBCellFAB&        a_primState,
        const Vector<VolIndex>& a_vofset,
        const int&              a_faceDir,
        const Side::LoHiSide&   a_sd,
        const Box&       a_box)
{
  //this holds velocity at time n
  const EBCellFAB& normalVel = *m_normalVelPtr;
  int nPrim = numPrimitives();
  for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real  vel = normalVel(vof, a_faceDir);

          for (int ivar = 0; ivar < nPrim; ivar++)
            {
              Real leftState, righState;
              if (a_sd == Side::Lo)
                {
                  leftState = a_exteState(vof, ivar);
                  righState = a_primState(vof, ivar);
                }
              else
                {
                  leftState = a_primState(vof, ivar);
                  righState = a_exteState(vof, ivar);
                }

              if (vel > TOL)
                {
                  a_coveredPrim(vof, ivar) = leftState;
                }
              else if (vel < -TOL)
                {
                  a_coveredPrim(vof, ivar) = righState;
                }
              else
                {
                  a_coveredPrim(vof, ivar) = 0.5*(righState+leftState);
                }
            }
        }
    }
}
/*****************************/
void EBPatchTransport::
getFlux(EBFluxFAB&       a_flux,
        const EBFluxFAB& a_prim)
{
  CH_assert(m_isAdvVelSet);
  int nFlux = numFluxes();
  CH_assert(nFlux == 1);

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {      
      const EBFaceFAB& advectVel = (*m_advVelPtr)[faceDir];
      a_flux[faceDir].setVal(0.);
      a_flux[faceDir].plus(a_prim[faceDir], 0, 0, nFlux);
      a_flux[faceDir].mult(advectVel, 0, 0, nFlux);
    }

}
/*****************************/
void EBPatchTransport::
updatePrim(EBCellFAB&              a_primMinu,
           EBCellFAB&              a_primPlus,
           const EBFaceFAB&        a_primFace,
           const BaseIVFAB<Real>&  a_coveredPrimMinu,
           const BaseIVFAB<Real>&  a_coveredPrimPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_faceDir,
           const Box&              a_box,
           const Real&             a_scale)
{
    //uses flux as prims
  //update the regular vofs.
  //first we have to copy the original state
  //so that the irregular stuff can be updated later.

  int numPrim = numPrimitives();
  //save state so that the regular can overwrite irregular cells
  Vector<Vector<Real> > cacheMinu;
  Vector<Vector<Real> > cachePlus;
  {
    CH_TIME("EBPatchAdvect::cache");
    cacheEBCF(cacheMinu, a_primMinu);
    cacheEBCF(cachePlus, a_primPlus);
  }

  const EBCellFAB& normalVel = *m_normalVelPtr;
  BaseFab<Real>& regPrimMinu = a_primMinu.getSingleValuedFAB();
  BaseFab<Real>& regPrimPlus = a_primPlus.getSingleValuedFAB();
  const BaseFab<Real>& regPrimFace      =    a_primFace.getSingleValuedFAB();

  {
    CH_TIME("EBPatchAdvect::fortran");
    const BaseFab<Real>& regNormVel = normalVel.getSingleValuedFAB();
    FORT_ADVECTUPDATE(CHF_BOX(a_box),
                      CHF_FRA(regPrimMinu),
                      CHF_FRA(regPrimPlus),
                      CHF_CONST_FRA(regPrimFace),
                      CHF_CONST_FRA(regNormVel),
                      CHF_CONST_INT(a_faceDir),
                      CHF_CONST_INT(numPrim),
                      CHF_CONST_REAL(a_scale),
                      CHF_BOX(a_box));
  }

  //put irregular cells back to their original state
  {
    CH_TIME("EBPatchAdvect::uncache");
    uncacheEBCF(a_primMinu ,cacheMinu);
    uncacheEBCF(a_primPlus ,cachePlus);
  }

  {
    CH_TIME("EBPatchAdvect::irregularLoop");
    for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        //the set can be bigger than the box for performance reasons.
        if (a_box.contains(vof.gridIndex()))
          {
            Real vel = normalVel(vof, a_faceDir);

            for (int ivar = 0; ivar < numPrim; ivar++)
              {

                Real primLo = 0;
                {
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, a_faceDir, Side::Lo);
                  if (facesLo.size() > 0)
                    {
                      for (int iface = 0; iface < facesLo.size(); iface++)
                        {
                          primLo += a_primFace(facesLo[iface], ivar);
                        }
                      primLo /= facesLo.size();
                    }
                  else
                    {
                      primLo = a_coveredPrimMinu(vof, ivar);
                    }
                }
                Real primHi = 0;
                {
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, a_faceDir, Side::Hi);
                  if (facesHi.size() > 0)
                    {
                      for (int iface = 0; iface < facesHi.size(); iface++)
                        {
                          primHi += a_primFace(facesHi[iface], ivar);
                        }
                      primHi /= facesHi.size();
                    }
                  else
                    {
                      primHi = a_coveredPrimPlus(vof, ivar);
                    }
                }
                Real primDiff = primHi - primLo;

                //state changed in regular update.
                //but we cached and uncached
                //scale holds 0.5*dt/dx in 2d.
                a_primMinu(vof, ivar) -=  a_scale*vel*primDiff;
                a_primPlus(vof, ivar) -=  a_scale*vel*primDiff;

                //debug set to inputs
                //a_primMinu(vof, ivar) = origPrimMinu;
                //a_primPlus(vof, ivar) = origPrimPlus;
                //end debug
              }
          }
      }
  }

  {
    CH_TIME("EBPatchAdvect::floors");
    floorPrimitives(a_primMinu, a_box);
    floorPrimitives(a_primPlus, a_box);
  } 
}
/*****************************/
EBPatchTransport::
EBPatchTransport():EBPatchGodunov()
{
  m_isNormalVelSet        = false;
  m_isAdvVelSet     = false;
  m_isMaxMinSet     = false;
  m_advVelPtr       = NULL;
  m_normalVelPtr    = NULL;
  m_useLimiting            = true;
  m_useFourthOrderSlopes   = false;
  m_useFlattening          = false;
}
/**************/
void
EBPatchTransport::
slope(EBCellFAB&       a_slopePrim,
      EBCellFAB&       a_slopeNLim,
      const EBCellFAB& a_primState,
      const EBCellFAB& a_flattening,
      const int&       a_dir,
      const Box&       a_box,
      bool a_doAgg)
{
  if (a_doAgg)
    {
      EBPatchGodunov::slope(a_slopePrim, a_slopeNLim, a_primState, a_flattening, a_dir, a_box, a_doAgg);
    }
  else
    {
      CH_assert(a_slopePrim.nComp() == numSlopes());
      CH_assert(a_primState.nComp() >= numSlopes());
      EBCellFAB deltaWL, deltaWR;

      {
        doSecondOrderSlopes(a_slopePrim,
                            deltaWL,
                            deltaWR,
                            a_slopeNLim,
                            a_primState,
                            a_dir,
                            a_box);
      }
    }
}
/*****************************/
void
EBPatchTransport::
doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                    EBCellFAB&       a_deltaWL,
                    EBCellFAB&       a_deltaWR,
                    EBCellFAB&       a_deltaWC,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const Box&       a_box)
{
  CH_assert(m_isSlopeSet);
  int numSlope = numSlopes();
  Box box1 = a_box;
  box1.grow(a_dir, 1);
  box1 &= m_domain;

  Box box2 = a_box;
  box2.grow(a_dir, 2);
  box2 &= m_domain;

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  if (usesFourthOrderSlopes())
    {
      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box2, m_domain, a_dir);
    }
  else
    {
      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box1, m_domain, a_dir);
    }

  {
    CH_TIME("defines_and_setvals");
    a_deltaWL.define(m_ebisBox, entireBox, numSlope);
    a_deltaWR.define(m_ebisBox, entireBox, numSlope);
    a_delta2W.setVal(0.);
    a_deltaWL.setVal(0.);
    a_deltaWR.setVal(0.);
    a_deltaWC.setVal(0.);
  }


  {
    CH_TIME("fortran");
    BaseFab<Real>& regDelta2W = a_delta2W.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWC = a_deltaWC.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWL = a_deltaWL.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWR = a_deltaWR.getSingleValuedFAB();
    const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();

    FORT_SECONDSLOPEDIFFS(CHF_FRA(regDeltaWC),
                          CHF_FRA(regDeltaWL),
                          CHF_FRA(regDeltaWR),
                          CHF_CONST_FRA(regPrimState),
                          CHF_CONST_INT(numSlope),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(loBox),
                          CHF_CONST_INT(hasLo),
                          CHF_BOX(hiBox),
                          CHF_CONST_INT(hasHi),
                          CHF_BOX(centerBox));

    //apply van leer limiter to regular cells.
    //limiting for irregular cells is more complicated
    //now that we are doing higher order limited
    //one-sided diffs
    regDelta2W.copy(regDeltaWC);
    if (m_useLimiting)
      {
        FORT_VLLIMITER(CHF_FRA(regDelta2W),
                       CHF_CONST_FRA(regDeltaWL),
                       CHF_CONST_FRA(regDeltaWR),
                       CHF_BOX(centerBox));
      }
  }

  {
    CH_TIME("irreg_slopes");
    //at irregular cells, if one side does not exist,
    //set all to one-sided diffs.  try to do higher order
    //limited diffs if possible.
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (entireBox.contains(vof.gridIndex()))
          {
            int ivar = 0;
            bool hasFacesLeft, hasFacesRigh;
            Real dql, dqr, dqc;
            bool verbose =false;
            pointGetSlopes(dql, dqr, dqc,
                           hasFacesLeft,
                           hasFacesRigh,
                           vof, a_primState, a_dir, ivar, verbose);
            Real dqSec = 0;
            if (!m_useLimiting)
              {
                dqSec = dqc;
              }
            else
              {
                if (hasFacesLeft && hasFacesRigh)
                  {
                    Real dqlim = dqc;

                    FORT_POINTVLLIMITER(CHF_REAL(dqlim),
                                        CHF_CONST_REAL(dql),
                                        CHF_CONST_REAL(dqr));
                    dqSec = dqlim;
                  }
                else if (!hasFacesLeft && !hasFacesRigh)
                  {
                    dqSec = 0.0;
                  }
                else if (hasFacesLeft && !hasFacesRigh)
                  {
                    if (dqc*dql > TOL)
                      {
                        Real rsign = 1.0;
                        if (dqc < -TOL)
                          {
                            rsign = -1.0;
                          }
                        dqSec = rsign*Min(Abs(dqc), Abs(dql));
                      }
                    else
                      {
                        dqSec = 0.0;
                      }

                  }
                else if (hasFacesRigh && !hasFacesLeft)
                  {
                    if (dqc*dqr > TOL)
                      {
                        Real rsign = 1.0;
                        if (dqc < -TOL)
                          {
                            rsign = -1.0;
                          }
                        dqSec = rsign*Min(Abs(dqc), Abs(dqr));
                      }
                    else
                      {
                        dqSec = 0.0;
                      }
                  }
                else
                  {
                    MayDay::Error("doh! missed a case");
                  }
              }//end if using limiting

            if (usesFlattening())
              {

                int pindex = pressureIndex();
                Real press = Max(a_primState(vof, pindex), Real(1.0e-10));
                Real dqp, dqlp, dqrp;
                verbose = false;
                pointGetSlopes(dqlp, dqrp, dqp,
                               hasFacesLeft,
                               hasFacesRigh,
                               vof, a_primState, a_dir, pindex,verbose);
                //try to detect pressure discontinutity
                dqp = Max(Abs(dqp), Abs(dqlp));
                dqp = Max(Abs(dqp), Abs(dqrp));
                if (Abs(dqp/press) > 0.1)
                  {
                    dqSec = 0.0;
                  }
              }
            a_delta2W(vof, ivar) = dqSec;
            a_deltaWL(vof, ivar) = dql;
            a_deltaWR(vof, ivar) = dqr;
            a_deltaWC(vof, ivar) = dqc;
          } //if vof is in this box
      }//end loop over vofs
  }
  {
    CH_TIME("bndry_slopes");
    //want to be able to call this in test codes
    if (m_isBCSet)
      {
        m_bc->setBndrySlopes(a_delta2W, a_primState, m_ebisBox, entireBox, a_dir);
      }
  }
}
/*****************************/
void
EBPatchTransport::
useLimiting(bool a_limiting)
{
  //sets flattening and fourth order slopes to false
  setSlopeParameters(false, false, a_limiting);
}
/*****************************/
EBPatchTransport::
~EBPatchTransport()
{
// already deleted
/*
  if (m_normalVelPtr != NULL) delete m_normalVelPtr;
  if (m_advVelPtr != NULL) delete m_advVelPtr;
  if (m_coveredAdvVelMinuPtr != NULL) delete m_coveredAdvVelMinuPtr;
  if (m_coveredAdvVelPlusPtr != NULL) delete m_coveredAdvVelPlusPtr;
*/
}
/*****************************/
Vector<string>
EBPatchTransport::
stateNames()
{
  Vector<string> retval(numConserved());
  for (int ivar = 0; ivar < numConserved(); ivar++)
    {
      char strname[100];
      sprintf(strname, "scalar%d",ivar);
      retval[ivar] = string(strname);
    }
  return retval;
}
/*****************************/
Vector<string>
EBPatchTransport::
primNames()
{
  Vector<string> retval(numConserved());
  for (int ivar = 0; ivar < numConserved(); ivar++)
    {
      char strname[100];
      sprintf(strname, "scalar%d",ivar);
      retval[ivar] = string(strname);
    }
  return retval;
}
/*****************************/
void
EBPatchTransport::
computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                   const EBCellFAB&  a_primState,
                   const EBCellFAB   a_slopePrim[SpaceDim],
                   const IntVectSet& a_irregIVS,
                   const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchTransport::computeEBIrregFlux");
  // loop over the cells with irregular faces and extrapolate
  // half-time variables to the faces using crude Taylor series
  // approximation.
  int numPrim = numPrimitives();
  const EBCellFAB& ccVel = *m_normalVelPtr;
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (m_validBox.contains(vof.gridIndex()))
        {
          IntVect iv = vof.gridIndex();
          RealVect bndryCentroid = m_ebisBox.bndryCentroid(vof);
          bndryCentroid *= m_dx;
          for (int comp=0; comp<numPrim; comp++)
            {
              Real traceVal = 0.0;
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  Real gradMult = bndryCentroid[dir] - 0.5*m_dt*ccVel(vof, dir);
                  traceVal -= gradMult*(a_slopePrim[dir])(vof, comp)/m_dx[dir];
                }
              traceVal += 0.5*m_dt*a_source(vof, comp);
              traceVal += a_primState(vof, comp);
              a_ebIrregFlux(vof, comp) = traceVal;
            } // end loop over components
        }
    } // end loop over vofs
}
/*****************************/
void EBPatchTransport::
setVelocities(const EBCellFAB& a_normalVel,
              const EBFluxFAB& a_advVel,
              const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelMinu,
              const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelPlus)
{
  setNormalVel(a_normalVel);
  setAdvVel(a_advVel, a_coveredAdvVelMinu, a_coveredAdvVelPlus);
}
/*****************************/
void EBPatchTransport::
setAdvVel(const EBFluxFAB& a_advVel,
          const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelMinu,
          const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelPlus)
{
  m_isAdvVelSet = true;
  m_advVelPtr = &a_advVel;
  m_coveredAdvVelMinuPtr = &a_coveredAdvVelMinu;
  m_coveredAdvVelPlusPtr = &a_coveredAdvVelPlus;
}
/*****************************/
void EBPatchTransport::
setNormalVel(const EBCellFAB& a_normalVel)
{
  m_isNormalVelSet = true;
  m_normalVelPtr = &a_normalVel;
}
/*****************************/
void EBPatchTransport::
averageVelToCC(EBCellFAB&                        a_normalVel,
               const EBFluxFAB&                  a_advectionVel,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloLo,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloHi,
               const Box&                        a_box) const
{
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceNormVel = a_advectionVel[faceDir];

      //treat every cell as regular
      BaseFab<Real>&       regCellVel = a_normalVel.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceNormVel.getSingleValuedFAB();
      FORT_AVEFACETOCELL(CHF_FRA(regCellVel),
                         CHF_CONST_FRA1(regFaceVel,0),
                         CHF_CONST_INT(faceDir),
                         CHF_BOX(a_box));
    }
  //correct way
  {
    CH_TIME("EBPatchAdvect::irregularLoop");
    for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        //the set can be bigger than the box for performance reasons.
        if (a_box.contains(vof.gridIndex()))
          {
            for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
              {
                Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, faceDir, Side::Lo);
                Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, faceDir, Side::Hi);
                Real velLo= 0.0;
                if (facesLo.size() > 0)
                  {
                    for (int iface = 0; iface < facesLo.size(); iface++)
                      {
                        velLo += a_advectionVel[faceDir](facesLo[iface], 0);
                      }
                    velLo /=  facesLo.size();
                  }
                else
                  {
                    const BaseIVFAB<Real>& coveredVelLo = *a_coveredVeloLo[faceDir];
                    velLo = coveredVelLo(vof, 0);
                  }
                Real velHi= 0.0;
                if (facesHi.size() > 0)
                  {
                    for (int iface = 0; iface < facesHi.size(); iface++)
                      {
                        velHi += a_advectionVel[faceDir](facesHi[iface], 0);
                      }
                    velHi /=  facesHi.size();
                  }
                else
                  {
                    const BaseIVFAB<Real>& coveredVelHi = *a_coveredVeloHi[faceDir];
                    velHi = coveredVelHi(vof, 0);
                  }

                Real velAve = 0.5*(velHi + velLo);
                a_normalVel(vof, faceDir) = velAve;
              }//end loop over directions
          }
      }//end loop over irregular vofs
  }
}  
/************/
void EBPatchTransport::
setCoveredConsVals(EBCellFAB& a_consState)
{
  a_consState.setInvalidData(0.0, 0);
}
/************/
/************/
void EBPatchTransport::
floorPrimitives(EBCellFAB&  a_primState,
                const Box&  a_box)
{
  if (m_isMaxMinSet)
    {
      CH_assert(numPrimitives()==1);
      BaseFab<Real>& primReg = a_primState.getSingleValuedFAB();

      FORT_FLOORPRIM(CHF_BOX(a_box),
                     CHF_REAL(m_minVal),
                     CHF_REAL(m_maxVal), 
                     CHF_FRA(primReg));

      for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
       {
         const VolIndex& vof = m_irregVoFs[ivof];
         if(a_box.contains(vof.gridIndex()))
           {
             a_primState(vof, 0) = Max(a_primState(vof, 0), m_minVal);
             a_primState(vof, 0) = Min(a_primState(vof, 0),  m_maxVal);
           }
       }
    }
}
/************/
void EBPatchTransport::
floorPrimitives(BaseIVFAB<Real>& a_primState,
                const IntVectSet& a_ivsIrreg)
{
  if (m_isMaxMinSet)
    {
      for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          //floors
          a_primState(vof, 0) = Min(a_primState(vof, 0),  m_maxVal);
          a_primState(vof, 0) = Max(a_primState(vof, 0),  m_minVal);
        }
    }
}
/************/
void EBPatchTransport::
floorConserved(EBCellFAB&  a_consState,
               const Box&  a_box)
{
  if (m_isMaxMinSet)
    {
      CH_assert(numConserved()==1);
      BaseFab<Real>& consReg = a_consState.getSingleValuedFAB();

      FORT_FLOORCONS(CHF_BOX(a_box),
                     CHF_REAL(m_minVal),
                     CHF_REAL(m_maxVal),
                     CHF_FRA(consReg));

      for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
       {
         const VolIndex& vof = m_irregVoFs[ivof];
         if(a_box.contains(vof.gridIndex()))
           {
             a_consState(vof, 0)  = Max(a_consState(vof, 0), m_minVal);
             a_consState(vof, 0) = Min(a_consState(vof, 0),  m_maxVal);
           }
        }
     }
}
/************/ 
void EBPatchTransport::
floorConserved(BaseIVFAB<Real>& a_consState,
               const IntVectSet& a_ivsIrreg)
{ 
  if (m_isMaxMinSet)
    { 
      for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        { 
          const VolIndex& vof = vofit();
          //floors 
          a_consState(vof, 0) = Min(a_consState(vof, 0),  m_maxVal);
          a_consState(vof, 0) = Max(a_consState(vof, 0),  m_minVal);
        } 
    }
} 
/************/
void EBPatchTransport::
primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                         EBCellFAB&          a_consState,
                         EBFluxFAB&          a_facePrim,
                         BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                         BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                         Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                         Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                         EBFluxFAB&          a_flux,
                         BaseIVFAB<Real>&    a_ebIrregFlux,
                         BaseIVFAB<Real>&    a_nonConservativeDivergence,
                         const EBCellFAB&    a_flattening,
                         const EBCellFAB&    a_source,
                         const Box&          a_box,
                         const IntVectSet&   a_ivsSmall,
                         const DataIndex&    a_dit,
                         bool                a_verbose)
{
  CH_assert(isDefined());
  int numCons = numConserved();

  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_consState.nComp() >= numConserved());
  setCoveredConsVals(a_consState);

  EBCellFAB slopePrim[SpaceDim];
  EBCellFAB slopeNLim[SpaceDim];

  a_facePrim.define(m_ebisBox, a_box, numPrimitives());

  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      IntVectSet ivsMinu = m_coveredSetsMinuG4[faceDir] & a_box;
      IntVectSet ivsPlus = m_coveredSetsPlusG4[faceDir] & a_box;
      a_coveredPrimMinu[faceDir].define(ivsMinu, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimPlus[faceDir].define(ivsPlus, m_ebisBox.getEBGraph(), numPrimitives());
    }

  extrapolatePrim(a_facePrim,
                  m_primState,
                  slopePrim,
                  slopeNLim,
                  a_coveredPrimMinu,
                  a_coveredPrimPlus,
                  a_coveredFaceMinu,
                  a_coveredFacePlus,
                  a_flattening,
                  a_consState,
                  a_source,
                  a_box,
                  a_dit,
                  a_verbose);

  getFlux(a_flux, a_facePrim); 

  if (usesArtificialViscosity())
    {
      MayDay::Error("EBPatchTransport::artificial viscosity is not expected");
    }

  // stable, non-conservative estimate of the solution update:
  advectiveDerivative(a_nonConsDivF, 
                      a_facePrim,
                      *m_advVelPtr,
                      a_coveredPrimMinu,
                      a_coveredPrimPlus,
                      *m_coveredAdvVelMinuPtr,
                      *m_coveredAdvVelPlusPtr,
                      a_coveredFaceMinu,
                      a_coveredFacePlus,
                      a_box);

  //copy nonconservative div f into sparse output thingy
  IntVectSet ivsIrreg = a_ivsSmall;
  for(VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for(int ivar = 0; ivar < numCons; ivar++)
        {
          a_nonConservativeDivergence(vof, ivar) = a_nonConsDivF(vof, ivar);
        }
    }

  //compute irregular boundary flux.  this includes
  //an artificial viscosity modification if appropriate
  computeEBIrregFlux( a_ebIrregFlux, m_primState,
                      slopeNLim, ivsIrreg, a_source);
}
/*******************/
void
EBPatchTransport::
extrapolatePrim(EBFluxFAB&                       a_facePrim,
                EBCellFAB&                       a_primState,
                EBCellFAB                        a_slopePrim[SpaceDim],
                EBCellFAB                        a_slopeNLim[SpaceDim],
                BaseIVFAB<Real>                  a_coveredPrimMinu[SpaceDim],
                BaseIVFAB<Real>                  a_coveredPrimPlus[SpaceDim],
                Vector<VolIndex>                 a_coveredFaceMinu[SpaceDim],
                Vector<VolIndex>                 a_coveredFacePlus[SpaceDim],
                const EBCellFAB&                 a_flattening,
                const EBCellFAB&                 a_consState,
                const EBCellFAB&                 a_source,
                const Box&                       a_box,
                const DataIndex&                 a_dit,
                bool                             a_verbose)
{
  #if CH_SPACEDIM==2

  extrapolatePrim2D(m_primMinu, m_primPlus,
                    m_primState, a_slopePrim, a_slopeNLim,
                    a_flattening, a_consState, a_source,
                    a_box, a_dit, a_verbose);


#elif CH_SPACEDIM==3

  extrapolatePrim3D(m_primMinu,   m_primPlus,
                    m_primState,  a_slopePrim, a_slopeNLim,
                    a_flattening, a_consState, a_source,
                    a_box, a_dit,a_verbose);

#else
  MayDay::Error("bogus SpaceDim");
#endif

  IntVectSet        coveredSetsPlus[SpaceDim];
  IntVectSet        coveredSetsMinu[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  //this keeps the fluxes from being calculated
  //on boundaries
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      bndryFaceBox[dir1] = a_box;
      bndryFaceBox[dir1] &= m_domain;
      bndryFaceBox[dir1].surroundingNodes(dir1);

      faceBox[dir1] = a_box;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVectSet irregIVSPlus, irregIVSMinu;
      computeCoveredFaces(a_coveredFacePlus[idir],
                          coveredSetsPlus[idir],
                          irregIVSPlus,
                          idir, Side::Hi, a_box);
      computeCoveredFaces(a_coveredFaceMinu[idir],
                          coveredSetsMinu[idir],
                          irregIVSMinu,
                          idir, Side::Lo, a_box);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_extendStatePlusG4[idir].setVal(0.);
      m_extendStateMinuG4[idir].setVal(0.);
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //extrapolate to covered faces using updated extrapolated state
      EBPatchGodunov::extrapToCoveredFaces(m_extendStateMinuG4[faceDir],
                                           m_primMinu[faceDir],
                                           m_primPlus[faceDir],
                                           m_primState,
                                           a_coveredFaceMinu[faceDir],
                                           faceDir, Side::Lo, a_box);

      EBPatchGodunov::extrapToCoveredFaces(m_extendStatePlusG4[faceDir],
                                           m_primMinu[faceDir],
                                           m_primPlus[faceDir],
                                           m_primState,
                                           a_coveredFacePlus[faceDir],
                                           faceDir, Side::Hi, a_box);


      riemann(a_facePrim[faceDir], m_primPlus[faceDir], m_primMinu[faceDir],
              faceDir, faceBox[faceDir]);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(a_facePrim, a_primState, m_primMinu[faceDir],
                   Side::Lo,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);
      m_bc->fluxBC(a_facePrim, a_primState, m_primPlus[faceDir],
                   Side::Hi,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);
      //solve riemann problem between extended state and
      //value in vof to get covered flux.
      riemann(a_coveredPrimMinu[faceDir],
              m_extendStateMinuG4[faceDir], m_primMinu[faceDir],
              a_coveredFaceMinu[faceDir], faceDir, Side::Lo, a_box);
      riemann(a_coveredPrimPlus[faceDir],
               m_extendStatePlusG4[faceDir], m_primPlus[faceDir],
              a_coveredFacePlus[faceDir], faceDir, Side::Hi, a_box);
    }
}
/****************/
void
EBPatchTransport::
advectiveDerivative(EBCellFAB&                      a_uDotDelRho,
                    const EBFluxFAB&                a_faceRho,
                    const EBFluxFAB&                a_faceVel,
                    const BaseIVFAB<Real>           a_coveredRhoLo[SpaceDim],
                    const BaseIVFAB<Real>           a_coveredRhoHi[SpaceDim],
                    const Vector<BaseIVFAB<Real>*>& a_coveredVelLo,
                    const Vector<BaseIVFAB<Real>*>& a_coveredVelHi,
                    const Vector<VolIndex>          a_coveredFaceLo[SpaceDim],
                    const Vector<VolIndex>          a_coveredFaceHi[SpaceDim],
                    const Box&                      a_box)
{
  int ncomp = a_faceRho.nComp();
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);

  CH_assert(a_coveredVelLo.size() == SpaceDim);
  CH_assert(a_coveredVelHi.size() == SpaceDim);

  //set udotdelrho to zero.  the fortran is additive
  a_uDotDelRho.setVal(0.);
  BaseFab<Real>& regUDotDelRho = a_uDotDelRho.getSingleValuedFAB();
  //compute udotdelu everywhere as regular.
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceRho = a_faceRho[faceDir];
      const EBFaceFAB& faceVel = a_faceVel[faceDir];

      const BaseFab<Real>& regFaceRho = faceRho.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceVel.getSingleValuedFAB();

      //this does udotdelvel += 0.5*(uhigh+ulow)*(velhigh-vellow)/dx (non-conservative)
      //this does udotdelrho +=   (velhigh*rhohigh-vellow*rholow)/dx (    conservative)
      //where high == ivec + half*e^facedir
      FORT_ADVECTIVEF(CHF_FRA(regUDotDelRho),
                      CHF_CONST_FRA(regFaceRho),
                      CHF_CONST_FRA1(regFaceVel, 0),
                      CHF_CONST_INT(faceDir),
                      CHF_CONST_INT(ncomp),
                      CHF_CONST_REAL(m_dx[faceDir]),
                      CHF_BOX(a_box),
                      CHF_INT(s_doingVel));
    }

  //update the irregular vofs
  for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      //the set can be bigger than the box for performance reasons.
      if (a_box.contains(vof.gridIndex()))
        {
          // if (s_verbose && vof.gridIndex()[0]==12 && vof.gridIndex()[1]==8){pout() << "nonconservative advection " << vof << endl;}
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              //udelrho was set in regular update.  we reset it
              // to zero and recalc.
              Real uDelRhoPt = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, idir, Side::Lo);
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, idir, Side::Hi);
                  Real rhoLo = 0.0;
                  Real velLo = 0.0;
                  if (facesLo.size() > 0)
                    {
                      for (int iface = 0; iface < facesLo.size(); iface++)
                        {
                          rhoLo += a_faceRho[idir](facesLo[iface], ivar);
                          velLo += a_faceVel[idir](facesLo[iface], 0);
                        }
                      velLo /=  facesLo.size();
                      rhoLo /=  facesLo.size();
                    }
                  else
                    {
                      const BaseIVFAB<Real>& coveredRhoLo = a_coveredRhoLo[idir];
                      const BaseIVFAB<Real>& coveredVelLo = *a_coveredVelLo[idir];
                      rhoLo = coveredRhoLo(vof, ivar);
                      velLo = coveredVelLo(vof, 0);
                    }
                  Real rhoHi = 0.0;
                  Real velHi = 0.0;
                  if (facesHi.size() > 0)
                    {
                      for (int iface = 0; iface < facesHi.size(); iface++)
                        {
                          rhoHi += a_faceRho[idir](facesHi[iface], ivar);
                          velHi += a_faceVel[idir](facesHi[iface], 0);
                        }
                      velHi /=  facesHi.size();
                      rhoHi /=  facesHi.size();
                    }
                  else
                    {
                      const BaseIVFAB<Real>& coveredRhoHi = a_coveredRhoHi[idir];
                      const BaseIVFAB<Real>& coveredVelHi = *a_coveredVelHi[idir];
                      rhoHi = coveredRhoHi(vof, ivar);
                      velHi = coveredVelHi(vof, 0);
                    }

                  Real velAve = 0.5*(velHi + velLo);
                  Real rhoDiff = rhoHi - rhoLo;
                  uDelRhoPt += velAve*rhoDiff/m_dx[idir];
                } //end loop over directions
              a_uDotDelRho(vof, ivar) = uDelRhoPt;
              // if (s_verbose && vof.gridIndex()[0]==12 && vof.gridIndex()[1]==8){pout() << "   uDotDelRho " << uDelRhoPt << endl;}
            }//end loop over variables
        } //if vof is in this box
    }//end loop over vofs 
}
/************/
#include "NamespaceFooter.H"
