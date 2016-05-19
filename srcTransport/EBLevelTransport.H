#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtg 1-29-2004

#ifndef _EBLEVELTRANSPORT_H_
#define _EBLEVELTRANSPORT_H_

#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "EBPWLFillPatch.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "EBPatchGodunovFactory.H"
#include "EBPatchTransport.H"
#include "EBPatchAdvect.H"
#include "NamespaceHeader.H"

///
/**
   Advection and diffusion on a level using
 */
class EBLevelTransport
{
public:
  /// Default constructor
  /**
     Object requires define(..) to be called before all other functions.
  */
  EBLevelTransport();

  /// Destructor
  /**
     destroys all objects created by define(..). Passed in
     data references of define(..) are left alone.
  */
  ~EBLevelTransport();

  EBLevelTransport(const DisjointBoxLayout&           a_thisDBL,
                   const DisjointBoxLayout&           a_coarDBL,
                   const EBISLayout&                  a_thisEBISL,
                   const EBISLayout&                  a_coarEBISL,
                   const ProblemDomain&               a_DProblem,
                   const int&                         a_nRefine,
                   const RealVect&                    a_dx,
                   const bool&                        a_hasCoarser,
                   const bool&                        a_hasFiner,
                   RefCountedPtr<EBPatchTransport>&   a_patchGodunov,
                   const EBIndexSpace*          const a_eb = Chombo_EBIS::instance());

  /// Actual constructor.
  /**
     For the coarsest level, an empty DisjointBoxLayout
     is passed in for coaserDisjointBoxLayout.
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const.
     (coding standards).  a_nRefine is the refinement ratio between
     this level and the next coarser level.  a_numGhosts is the number
     of ghost cells in each direction.
  */
  void define(const DisjointBoxLayout&           a_thisDBL,
              const DisjointBoxLayout&           a_coarDBL,
              const EBISLayout&                  a_thisEBISL,
              const EBISLayout&                  a_coarEBISL,
              const ProblemDomain&               a_DProblem,
              const int&                         a_nRefine,
              const RealVect&                    a_dx,
              const bool&                        a_hasCoarser,
              const bool&                        a_hasFiner,
              RefCountedPtr<EBPatchTransport>&   a_patchGodunov,
              const EBIndexSpace*                const a_eb = Chombo_EBIS::instance());

  bool isDefined() const;

  void floorPrimitives(LevelData<EBCellFAB>&         a_primState,
                       Real a_time, Real a_dt);

  void floorConserved(LevelData<EBCellFAB>&         a_consState,
                       Real a_time, Real a_dt);

  /// outputs divergeF.  a_cons gets its ghost cells filled
  /**
     For the coarsest level, UcoarseOld and UCoarseNew are empty
     LevelData<EBCellFAB> objects.  Also, either UCoarseOld or
     UCoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  U must be defined on the same
     DisjointBoxLayouts as were used in define(...)
     Coarse flux register is flux register with the coarse level.
     Fine flux register is the flux register with the fine level.
     To the finer level FR, this level is the coarse level.
     To the coarser level FR, this level is the fine level.
     outputs divergeF.  a_cons gets its ghost cells filled
  */
  void divergeF(LevelData<EBCellFAB>&         a_divergeF,
                LevelData<BaseIVFAB<Real> >&  a_massDiff,
                EBFluxRegister&               a_fineFluxRegister,
                EBFluxRegister&               a_coarFluxRegister,
                LevelData<EBCellFAB>&         a_consState,   //not really changed
                LevelData<EBCellFAB>&         a_normalVel,
                const LevelData<EBFluxFAB>&   a_advVel,
                const LayoutData< Vector <BaseIVFAB<Real> * > >& a_coveredAdvVelMinu,
                const LayoutData< Vector <BaseIVFAB<Real> * > >& a_coveredAdvVelPlus,
                const LevelData<EBCellFAB>&   a_source,
                const LevelData<EBCellFAB>&   a_consStateCoarseOld,
                const LevelData<EBCellFAB>&   a_consStateCoarseNew,
                const LevelData<EBCellFAB>&   a_normalVelCoarseOld,
                const LevelData<EBCellFAB>&   a_normalVelCoarseNew,
                const Real&                   a_time,
                const Real&                   a_coarTimeOld,
                const Real&                   a_coarTimeNew,
                const Real&                   a_dt);

  void computeNormalVel(LevelData<EBCellFAB>&                   a_normalVel,
                 const LevelData<EBFluxFAB>&                    a_advectionVel,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloLo,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloHi); 

protected:

  void fillCons(LevelData<EBCellFAB>&         a_consState,
                LevelData<EBCellFAB>&         a_normalVel,
                const LevelData<EBCellFAB>&   a_consStateCoarseOld,
                const LevelData<EBCellFAB>&   a_consStateCoarseNew,
                const LevelData<EBCellFAB>&   a_normalVelCoarseOld,
                const LevelData<EBCellFAB>&   a_normalVelCoarseNew,
                const Real&                   a_time,
                const Real&                   a_coarTimeOld,
                const Real&                   a_coarTimeNew);

  void computeFlattening(LevelData<EBCellFAB>&         a_cons,
                         Real a_time, Real a_dt);

  void doRegularUpdate(LevelData<EBCellFAB>&         a_divergeF,
                       LevelData<EBCellFAB>&         a_cons,
                       const LevelData<EBCellFAB>&   a_normalVel,
                       const LevelData<EBFluxFAB>&   a_advVel,
                       const LayoutData< Vector <BaseIVFAB<Real> * > >& a_coveredAdvVelMinu,
                       const LayoutData< Vector <BaseIVFAB<Real> * > >& a_coveredAdvVelPlus,
                       EBFluxRegister&               a_fineFluxRegister,
                       EBFluxRegister&               a_coarFluxRegister,
                       const LevelData<EBCellFAB>&   a_source,
                       Real a_time, Real a_dt);

  void doIrregularUpdate(LevelData<EBCellFAB>&         a_divergeF,
                         LevelData<EBCellFAB>&         a_cons,
                         EBFluxRegister&               a_fineFluxRegister,
                         EBFluxRegister&               a_coarFluxRegister,
                         LevelData<BaseIVFAB<Real> >&  a_massDiff,
                         Real a_time, Real a_dt);

  bool                       m_hasCoarser;
  bool                       m_hasFiner;
  bool                       m_isDefined;
  RealVect                   m_dx;
  int                        m_refRatCrse;
  int                        m_nVar;
  int                        m_nPrim;
  int                        m_nCons;
  int                        m_nFlux;
  int                        m_nGhost;

  //these are not grown by one.
  LayoutData<IntVectSet> m_irregSetsSmall;

  //these are grown by one in the directions != idir
  //gets filled in EBArith::defineFluxInterpolant
  LayoutData<IntVectSet> m_irregSetsGrown[SpaceDim];

  LevelData<BaseIFFAB<Real> > m_fluxInterpolants[SpaceDim];
  LevelData<BaseIVFAB<Real> >    m_nonConsDivergence;
  LevelData<BaseIVFAB<Real> >    m_ebIrregFaceFlux;
  LayoutData<IntVectSet>     m_cfIVS;
  DisjointBoxLayout          m_thisGrids;
  DisjointBoxLayout          m_coarGrids;
  EBISLayout                 m_thisEBISL;
  EBISLayout                 m_coarEBISL;
  EBPWLFillPatch             m_fillPatch;
  EBPWLFillPatch             m_fillPatchVel;
  ProblemDomain              m_domain;
//  LayoutData<EBPatchTransport*> m_ebPatchTransport;
  RefCountedPtr<EBPatchTransport>  m_ebPatchGodunov;
//  RefCountedPtr<EBPatchAdvect>  m_ebPatchGodunov;
  LevelData<EBCellFAB>          m_flattening; 

private:
  //disallowed for all the usual reasons
  void operator=(const EBLevelTransport& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBLevelTransport(const EBLevelTransport& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
