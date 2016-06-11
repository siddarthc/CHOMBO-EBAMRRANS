#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "kEpsilon.H"

#include "NamespaceHeader.H"
/**********/
kEpsilon::
kEpsilon():EBPatchRANSModel()
{
  setNumEquations();

// A good place to hardwire I think
// But this may need changing in derivarive's constructors. Beware!
  m_nPrim = m_nEqn;
  m_nCons = m_nEqn;
  m_nFlux = m_nEqn;
  m_nSlope = m_nEqn;
  m_kEqnIndex = 0;
  m_epsEqnIndex = 1;
}
/**********/
kEpsilon::
~kEpsilon()
{

}
/**********/
void kEpsilon::setNumEquations()
{
  m_nEqn = 2;
}
/**********/
void kEpsilon::floorPrimitives(EBCellFAB& a_primState,
                               const Box& a_box)
{
 // feeling lazy; just using grandparents stuff
  if (m_params.m_limitK)
  {
    EBCellFAB tempFAB(m_ebisBox,a_box,1);
    Interval kInterv(m_kEqnIndex,m_kEqnIndex);
    Interval tempInterv(0,0);
    tempFAB.copy(a_box,tempInterv,a_box,a_primState,kInterv);
    EBPatchTransport::setMaxMin(m_params.m_maxK,m_params.m_minK);
    EBPatchTransport::floorPrimitives(tempFAB,a_box);
    a_primState.copy(a_box,kInterv,a_box,tempFAB,tempInterv);
  }
//  m_isMaxMinSet = false;
  if (m_params.m_limitEpsilon)
  {
    EBCellFAB tempFAB(m_ebisBox,a_box,1);
    Interval epsInterv(m_epsEqnIndex,m_epsEqnIndex);
    Interval tempInterv(0,0);
    tempFAB.copy(a_box,tempInterv,a_box,a_primState,epsInterv);
    EBPatchTransport::setMaxMin(m_params.m_maxEpsilon,m_params.m_minEpsilon);
    EBPatchTransport::floorPrimitives(tempFAB,a_box);
    a_primState.copy(a_box,epsInterv,a_box,tempFAB,tempInterv);
  }  
}
/**********/
void kEpsilon::floorPrimitives(BaseIVFAB<Real>&  a_primState,
                                const IntVectSet& a_ivsIrreg)
{
  if (m_params.m_limitK) 
  {
    for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          a_primState(vof,m_kEqnIndex) = Min(a_primState(vof,m_kEqnIndex), m_params.m_minK);
          a_primState(vof,m_kEqnIndex) = Max(a_primState(vof,m_kEqnIndex), m_params.m_maxK);
        }
  }
//  m_isMaxMinSet = false;
  if (m_params.m_limitEpsilon)
  {
    for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          a_primState(vof,m_epsEqnIndex) = Min(a_primState(vof,m_epsEqnIndex), m_params.m_minEpsilon);
          a_primState(vof,m_epsEqnIndex) = Max(a_primState(vof,m_epsEqnIndex), m_params.m_maxEpsilon); 
        }
  }
}
/**********/
void kEpsilon::floorConserved(EBCellFAB&  a_consState,
                              const Box&  a_box)
{
 // feeling lazy; just using grandparents stuff
  if (m_params.m_limitK)
  {
    EBCellFAB tempFAB(m_ebisBox,a_box,1);
    Interval kInterv(m_kEqnIndex,m_kEqnIndex);
    Interval tempInterv(0,0);
    tempFAB.copy(a_box,tempInterv,a_box,a_consState,kInterv);
    EBPatchTransport::setMaxMin(m_params.m_maxK,m_params.m_minK);
    EBPatchTransport::floorConserved(tempFAB,a_box);
    a_consState.copy(a_box,kInterv,a_box,tempFAB,tempInterv);
  }
//  m_isMaxMinSet = false;
  if (m_params.m_limitEpsilon)
  {
    EBCellFAB tempFAB(m_ebisBox,a_box,1);
    Interval epsInterv(m_epsEqnIndex,m_epsEqnIndex);
    Interval tempInterv(0,0);
    tempFAB.copy(a_box,tempInterv,a_box,a_consState,epsInterv);
    EBPatchTransport::setMaxMin(m_params.m_maxEpsilon,m_params.m_minEpsilon);
    EBPatchTransport::floorConserved(tempFAB,a_box);
    a_consState.copy(a_box,epsInterv,a_box,tempFAB,tempInterv);
  }
}
/**********/
void kEpsilon::floorConserved(BaseIVFAB<Real>&   a_consState,
                              const IntVectSet&  a_ivsIrreg)
{
  if (m_params.m_limitK)
  {
    for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          a_consState(vof,m_kEqnIndex) = Min(a_consState(vof,m_kEqnIndex), m_params.m_minK);
          a_consState(vof,m_kEqnIndex) = Max(a_consState(vof,m_kEqnIndex), m_params.m_maxK);
        }
  }
//  m_isMaxMinSet = false;
  if (m_params.m_limitEpsilon)
  {
    for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          a_consState(vof,m_epsEqnIndex) = Min(a_consState(vof,m_epsEqnIndex), m_params.m_minEpsilon);
          a_consState(vof,m_epsEqnIndex) = Max(a_consState(vof,m_epsEqnIndex), m_params.m_maxEpsilon);
        }
  }
}
/**********/
Vector<string>
kEpsilon::stateNames()
{
  Vector<string> retVal(m_nEqn);

  retVal[m_kEqnIndex] = "k";
  retVal[m_epsEqnIndex] = "epsilon";
  
  return retVal;  
}
/**********/
Vector<string>
kEpsilon::primNames()
{
  Vector<string> retVal(m_nEqn);
  
  retVal[m_kEqnIndex] = "k";
  retVal[m_epsEqnIndex] = "epsilon";
  
  return retVal;   
}
/***********/
void kEpsilon::
computeNetSource(EBCellFAB&       a_netSource,
                 const EBCellFAB& a_state,
                 const Box&       a_box)
{

}
/**********/
Real kEpsilon::
stencilWeightEqnIndex()
{
  // reset stencil weights based on kEqn
  return m_kEqnIndex;
}
/*********/
void kEpsilon::
getDiffusionCoefficients(EBFluxFAB&       a_diffCoeff,
                         const EBFluxFAB& a_stateFace,
                         const Box&       a_box)
{

}
/**********/
void kEpsilon::
getDiffusionCoefficients(BaseIVFAB<Real>& a_diffCoeffIrreg,
                         const EBCellFAB& a_stateCell,
                         const IntVectSet& a_ivs)
{

}
/*********/
void kEpsilon::
fillDriverDiffusionCoefficients(EBFluxFAB&       a_diffCoeff,
                         const EBFluxFAB& a_stateFace,
                         const Box&       a_box)
{

}
/**********/
void kEpsilon::
fillDriverDiffusionCoefficients(BaseIVFAB<Real>& a_diffCoeffIrreg,
                         const EBCellFAB& a_stateCell,
                         const IntVectSet& a_ivs)
{

}
/**********/
#include "NamespaceFooter.H"
