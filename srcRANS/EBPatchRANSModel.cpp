#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchRANSModel.H"
#include "EBPatchRANSModelF_F.H"
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
/***************/
EBPatchRANSModel::
EBPatchRANSModel():EBPatchTransport()
{
//  setNumEquations(); it'd be nice if i can do this
//                     but it's a recipe for disaster 
//                     better way is to call this from factory
}
/***************/
EBPatchRANSModel::
~EBPatchRANSModel()
{
// already deleted
/*
  if (m_normalVelPtr != NULL) delete m_normalVelPtr;
  if (m_advVelPtr != NULL) delete m_advVelPtr;
  if (m_coveredAdvVelMinuPtr != NULL) delete m_coveredAdvVelMinuPtr;
  if (m_coveredAdvVelPlusPtr != NULL) delete m_coveredAdvVelPlusPtr;
*/  
}
/******/
void EBPatchRANSModel::
normalPred(EBCellFAB&       a_primLo,
           EBCellFAB&       a_primHi,
           const EBCellFAB& a_primState,
           const EBCellFAB& a_slopePrim,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_assert(m_isNormalVelSet);
  CH_assert(m_normalVelPtr != NULL);
  const EBCellFAB& velcc = *m_normalVelPtr;

    /**/
  BaseFab<Real>& regPrimLo    = a_primLo.getSingleValuedFAB();
  BaseFab<Real>& regPrimHi    = a_primHi.getSingleValuedFAB();
  const BaseFab<Real>& regVel =  velcc.getSingleValuedFAB();
  const BaseFab<Real>& regPrim =  a_primState.getSingleValuedFAB();
  const BaseFab<Real>& regSlope = a_slopePrim.getSingleValuedFAB();

  //saved multiple box iterations in fortran!
  FORT_PREDADVECTRANS(CHF_BOX(a_box),
                      CHF_CONST_FRA(regPrim),
                      CHF_CONST_FRA(regSlope),
                      CHF_CONST_FRA(regVel),
                      CHF_FRA(regPrimLo),
                      CHF_FRA(regPrimHi),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_REAL(a_dtbydx),
                      CHF_CONST_INT(m_nEqn));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> primlo(m_nEqn), primhi(m_nEqn);
      Vector<Real> primit(m_nEqn), pslope(m_nEqn);
      RealVect veloc;

      for (int ivar = 0; ivar < m_nEqn; ivar++)
        {
          primit[ivar] = a_primState(vof, ivar);
          pslope[ivar] = a_slopePrim(vof, ivar);
        }

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = velcc(vof, idir);
        }

      FORT_POINTPREDADVECTRANS(CHF_VR(primit),
                               CHF_VR(primlo),
                               CHF_VR(primhi),
                               CHF_VR(pslope),
                               CHF_REALVECT(veloc),
                               CHF_CONST_INT(a_dir),
                               CHF_CONST_REAL(a_dtbydx),
                               CHF_CONST_INT(m_nEqn));

      for (int ivar=0; ivar < m_nEqn; ivar++)
        {
          a_primLo(vof, ivar) = primlo[ivar];
          a_primHi(vof, ivar) = primhi[ivar];
        } 
    }  
}
/****************/
void EBPatchRANSModel::
getFlux(EBFluxFAB&       a_flux,
        const EBFluxFAB& a_prim)
{
  // random checks:
  CH_assert(m_isAdvVelSet);
  int nFlux = numFluxes();
  CH_assert(nFlux = m_nEqn);

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& advectVel = (*m_advVelPtr)[faceDir];
      a_flux[faceDir].setVal(0.);
      a_flux[faceDir].plus(a_prim[faceDir], 0, 0, nFlux);
      a_flux[faceDir].mult(advectVel, 0, 0, nFlux);
    }
}
/****************/
#include "NamespaceFooter.H"
