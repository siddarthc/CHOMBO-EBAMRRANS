#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBLevelRANS.H"
#include "EBPatchTransport.H"

#include "NamespaceHeader.H"

/***********/
EBLevelRANS::
EBLevelRANS() : EBLevelTransport()
{

}
/**********/
EBLevelRANS::
~EBLevelRANS()
{
}
/**********/
EBLevelRANS::
EBLevelRANS(const DisjointBoxLayout&           a_thisDBL,
            const DisjointBoxLayout&           a_coarDBL,
            const EBISLayout&                  a_thisEBISL,
            const EBISLayout&                  a_coarEBISL,
            const ProblemDomain&               a_domain,
            const int&                         a_nRefine,
            const RealVect&                    a_dx,
            const bool&                        a_hasCoarser,
            const bool&                        a_hasFiner,
            RefCountedPtr<EBPatchRANSModel>&   a_patchRANS,
            const EBIndexSpace*          const a_eb)
{
  CH_TIME("EBLevelRANS::EBLevelRANS");
  m_isDefined = false;
  define(a_thisDBL,
         a_coarDBL,
         a_thisEBISL,
         a_coarEBISL,
         a_domain,
         a_nRefine,
         a_dx,
         a_hasCoarser,
         a_hasFiner,
         a_patchRANS,
         a_eb);
}
/**********/
void 
EBLevelRANS::
define(const DisjointBoxLayout&           a_thisDBL,
       const DisjointBoxLayout&           a_coarDBL,
       const EBISLayout&                  a_thisEBISL,
       const EBISLayout&                  a_coarEBISL,
       const ProblemDomain&               a_DProblem,
       const int&                         a_nRefine,
       const RealVect&                    a_dx,
       const bool&                        a_hasCoarser,
       const bool&                        a_hasFiner,
       RefCountedPtr<EBPatchRANSModel>&   a_patchRANS,
       const EBIndexSpace*                const a_eb)
{
  m_ebPatchModel = a_patchRANS; 

// dangerous casting but dynamic cast complains about RefCountedPtr
// but I'm sure m_ebPatchModel is a EBPatchTransport derivative 
  RefCountedPtr<EBPatchTransport> patchTransport = static_cast<RefCountedPtr<EBPatchTransport> >(m_ebPatchModel);

  EBLevelTransport::define(a_thisDBL,a_coarDBL,a_thisEBISL,a_coarEBISL,a_DProblem,a_nRefine,a_dx,a_hasCoarser,a_hasFiner,patchTransport,a_eb);
}
/**********/
void 
EBLevelRANS::
computeNetSource(LevelData<EBCellFAB>&       a_netSource,
                 const LevelData<EBCellFAB>& a_state,
                 Real a_time, Real a_dt)
{
    for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& netSource = a_netSource[dit()];
      const EBCellFAB& primState = a_state[dit()];
      const IntVectSet& cfivs = m_cfIVS[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      const Box& cellBox = m_thisGrids.get(dit());
      m_ebPatchModel->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
      m_ebPatchModel->computeNetSource(netSource, primState, cellBox);
    }
}
/*********/
void 
EBLevelRANS::
getDiffusionCoefficients(LevelData<EBFluxFAB>&        a_diffCoeff,
                         LevelData<BaseIVFAB<Real> >& a_diffCoeffIrreg,
                         const LevelData<EBCellFAB>&  a_stateCell,
                         const LevelData<EBFluxFAB>&  a_stateFace,
                         Real a_time, Real a_dt)
{
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
   {
     EBFluxFAB& diffCoeff = a_diffCoeff[dit()];
     BaseIVFAB<Real>& diffCoeffIrreg = a_diffCoeffIrreg[dit()];
     const EBCellFAB& stateCell = a_stateCell[dit()];
     const EBFluxFAB& stateFace = a_stateFace[dit()];
     
     const IntVectSet& cfivs = m_cfIVS[dit()];
     const EBISBox& ebisBox = m_thisEBISL[dit()];
     const Box& cellBox = m_thisGrids.get(dit());
     m_ebPatchModel->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
     m_ebPatchModel->getDiffusionCoefficients(diffCoeff, stateFace, cellBox);

     const IntVectSet& ivs = ebisBox.getIrregIVS(cellBox);
     m_ebPatchModel->getDiffusionCoefficients(diffCoeffIrreg, stateCell, ivs);
   }  
}
/**********/
void EBLevelRANS::
fillDriverDiffusionCoefficients(RefCountedPtr<LevelData<EBFluxFAB> >&        a_diffCoeff,
                         RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_diffCoeffIrreg,
                         const LevelData<EBCellFAB>&  a_stateCell,
                         const LevelData<EBFluxFAB>&  a_stateFace,
                         Real a_time, Real a_dt)
{
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
   {
     EBFluxFAB& diffCoeff = (*a_diffCoeff)[dit()];
     BaseIVFAB<Real>& diffCoeffIrreg = (*a_diffCoeffIrreg)[dit()];
     const EBCellFAB& stateCell = a_stateCell[dit()];
     const EBFluxFAB& stateFace = a_stateFace[dit()];

     const IntVectSet& cfivs = m_cfIVS[dit()];
     const EBISBox& ebisBox = m_thisEBISL[dit()];
     const Box& cellBox = m_thisGrids.get(dit());
     m_ebPatchModel->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
     m_ebPatchModel->fillDriverDiffusionCoefficients(diffCoeff, stateFace, cellBox);
 
     const IntVectSet& ivs = ebisBox.getIrregIVS(cellBox);
     m_ebPatchModel->fillDriverDiffusionCoefficients(diffCoeffIrreg, stateCell, ivs);
   }
}
/**********/
#include "NamespaceFooter.H"
