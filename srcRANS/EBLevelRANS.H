#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELRANS_H_
#define _EBLEVELRANS_H_

#include "EBLevelTransport.H"
#include "EBPatchRANSModelFactory.H"
#include "EBPatchRANSModel.H"

#include "NamespaceHeader.H"

class EBLevelRANS : public EBLevelTransport
{
public:

  EBLevelRANS();

  ~EBLevelRANS();

  EBLevelRANS(const DisjointBoxLayout&           a_thisDBL,
              const DisjointBoxLayout&           a_coarDBL,
              const EBISLayout&                  a_thisEBISL,
              const EBISLayout&                  a_coarEBISL,
              const ProblemDomain&               a_DProblem,
              const int&                         a_nRefine,
              const RealVect&                    a_dx,
              const bool&                        a_hasCoarser,
              const bool&                        a_hasFiner,
              RefCountedPtr<EBPatchRANSModel>&   a_patchRANS,
              const EBIndexSpace*          const a_eb = Chombo_EBIS::instance());

  virtual void define(const DisjointBoxLayout&           a_thisDBL,
                      const DisjointBoxLayout&           a_coarDBL,
                      const EBISLayout&                  a_thisEBISL,
                      const EBISLayout&                  a_coarEBISL,
                      const ProblemDomain&               a_DProblem,
                      const int&                         a_nRefine,
                      const RealVect&                    a_dx,
                      const bool&                        a_hasCoarser,
                      const bool&                        a_hasFiner,
                      RefCountedPtr<EBPatchRANSModel>&   a_patchRANS,
                     const EBIndexSpace*                const a_eb = Chombo_EBIS::instance());  

  /// compute production - dissipation for each equation in the model
  virtual void computeNetSource(LevelData<EBCellFAB>&       a_netSource,
                                const LevelData<EBCellFAB>& a_state,
                                Real a_time, Real a_dt);

  ///set the equation diffusion coefficients
  virtual void getDiffusionCoefficients(LevelData<EBFluxFAB>&        a_diffCoeff, 
                                        LevelData<BaseIVFAB<Real> >& a_diffCoeffIrreg,
                                        const LevelData<EBCellFAB>&  a_stateCell, 
                                        const LevelData<EBFluxFAB>&  a_stateFace,
                                        Real a_time, Real a_dt);

protected:

  RefCountedPtr<EBPatchRANSModel> m_ebPatchModel;
};

#include "NamespaceFooter.H"
#endif
