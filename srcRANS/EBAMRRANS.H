#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRRANS_H_
#define _EBAMRRANS_H_

#include "EBAMRLevel.H"
#include "EBAMRRANSParams.H"
#include "EBPatchRANSModelFactory.H"

#include "EBCellFAB.H"
#include "EBLevelTGA.H"
#include "BaseIVFAB.H"
#include "LevelData.H"
#include "EBLevelRedist.H"
#include "EBCoarToCoarRedist.H"
#include "EBPatchGodunov.H"
#include "EBPatchGodunovFactory.H"
#include "EBBackwardEuler.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarseAverage.H"
#include "EBPWLFineInterp.H"
#include "EBFluxRegister.H"
#include "Box.H"
#include "BiCGStabSolver.H"
#include "EBSimpleSolver.H"
#include "IntVectSet.H"
#include "AMRMultiGrid.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "EBQuadCFInterp.H"
#include "BaseDomainBC.H"
#include "BaseEBBC.H"

#include "NamespaceHeader.H"

// RANS solver 
// driven by external AMR or EBAMR class
// I could derive this from EBAMRTransport to save time and effort but I want RANS to have its own thing

// For now this is only for incompressible RANS
//   but can be extended to compressible RANS
//     with modifications: variable aco, add routine getPrimState and get coefficients from primState
//       or may be there's a way to get away by making these corrections in the level and patch stuff

class EBAMRRANS : public EBAMRLevel
{
public:

  virtual ~EBAMRRANS();

  EBAMRRANS(const EBAMRRANSParams& a_params,
            const RefCountedPtr<EBPatchRANSModelFactory>& a_modelFactory,
            bool  a_externalDriver);

  virtual void define(EBAMRLevel*            a_coarser_level_ptr,
                      const ProblemDomain& a_problem_domain,
                      int                  a_level,
                      int                  a_ref_ratio);

  virtual void assignDx(RealVect a_dx);

  virtual Real advance();

  virtual void postTimeStep();

  virtual void tagCells(IntVectSet& a_tags);

  virtual void tagCellsInit(IntVectSet& a_tags);

  virtual void regrid(const Vector<Box>& a_new_grids);

  virtual void preRegrid(int                         a_base_level,
                         const Vector<Vector<Box> >& a_new_grids);

  virtual void postRegrid(int a_base_level);

  void defineSolvers();

  static bool s_solversDefined;

#ifdef CH_USE_HDF5
  // ---------------------------------------------------
  // HDF5 routines for plotting and writing checkpoints.
  // ---------------------------------------------------

  virtual void writePlotHeaderOld    (HDF5Handle& a_handle) const;
  virtual void writePlotLevelOld     (HDF5Handle& a_handle) const;
  virtual void writePlotHeader       (HDF5Handle& a_handle) const;
  virtual void writePlotLevel        (HDF5Handle& a_handle) const;
  virtual void writeCheckpointHeader (HDF5Handle& a_handle) const;
  virtual void writeCheckpointLevel  (HDF5Handle& a_handle) const;
  virtual void readCheckpointHeader  (HDF5Handle& a_handle);
  virtual void readCheckpointLevel   (HDF5Handle& a_handle);
#endif

protected:

  RefCountedPtr<EBPatchRANSModel> m_ebPatchModel;  

  void clearSolvers();
 
  void fluxDivergence(  LevelData<EBCellFAB>& a_divergeF);

  void hyperbolicSource(LevelData<EBCellFAB>& a_source);

  //this includes news so need to call delete
  void
  getCoarserState(LevelData<EBCellFAB>* & a_stateCoar);

  //if doNormalization = true, averages with neighboring cells
  //otherwise returns kappa * source
  void explicitHyperbolicSource(LevelData<EBCellFAB>&       a_hypSource,
                                const LevelData<EBCellFAB>& a_state,
                                bool a_doNormalization ); 

  /// set  output to volfrac*((del dot (kappa grad S))
  void kappaDiffusiveSource(LevelData<EBCellFAB>& a_kappaDiffSource,
                            const LevelData<EBCellFAB>& a_state,
                            const LevelData<EBCellFAB>* a_stateCoar); 

  /// a_kappaNetSource = volfrac*(production - dissipation)
  void kappaNetSource(LevelData<EBCellFAB>& a_kappaNetSource,
                      const LevelData<EBCellFAB>& a_state,
                      const LevelData<EBCellFAB>* a_stateCoar);

  void coarseFineIncrement();

  // inviscid advance
  void explicitAdvance(const LevelData<EBCellFAB>& a_divergeF);

  //this one is explicit and also does reweighting of the 
  //redistribution object in the case of mass weighting
  void hyperbolicRedistribution(LevelData<EBCellFAB>& a_state);

  void  getUStar(LevelData<EBCellFAB>      & a_UStar,
                 const LevelData<EBCellFAB>& a_UN,
                 const LevelData<EBCellFAB>& a_divergef);

  void addNetSource(LevelData<EBCellFAB>& a_netSource,
                    LevelData<EBCellFAB>& a_UStar); 

  void updateStateByNetSourceAndRedistribute(LevelData<EBCellFAB>& a_netSource,
                                             LevelData<EBCellFAB>& a_kappaConsNetSource,
                                             LevelData<EBCellFAB>& a_nonConsNetSource,
                                             LevelData<EBCellFAB>& a_UStar);

  void getDivDGradS(LevelData<EBCellFAB>& a_divDGradS,
                    LevelData<EBCellFAB>& a_Ustar);

  //put ustar into state
  void finalAdvance(LevelData<EBCellFAB>& a_Ustar);
 
  void postTimeStepRefluxRedistDance();

  void resetWeights();

  void refluxRedistInteraction();

  void coarseFineRedistribution(const Interval& a_interv);

  void explicitReflux(const Interval& a_interv);

  // This performs the redistribution of a quantity (mass, momentum, or energy) in m_stateNew. 
  // The quantity is determined by a_interv. If a_putIntoState is set, the redistribution is 
  // done EXPLICITLY, and the redistributed quantity is placed directly into m_stateNew.
  // Otherwise, the quantity is accumulated into m_redisRHS for the later implicit 
  // redistribution.
  void refluxRHS();

  void implicitReflux();

  int getFinestLevel();

  //(alpha I - dt Ls) delta = dt*Dr(Fs)
  void getRefluxDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                       Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeS,
                       int baseLev, int finestLev, Real baseDt);

  //scalar += dt*rho*Cv*(deltaT)
  void incrScalarByDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                          int baseLev, int finestLev);
  void incrScalarByDeltaS(LevelData<EBCellFAB>& a_state,
                          LevelData<EBCellFAB>& a_deltaScalar);

  void defineFactories(bool a_atHalfTime);

  void getHalfState(LevelData<EBCellFAB>& a_stateInt);
  
  void fillCoefficients(const LevelData<EBCellFAB>& a_state);

  void setBCs();

  // no. of equations in the model
  int m_nEqn;

  bool m_externalDriver;
  EBAMRRANSParams m_params;
  RefCountedPtr<EBPatchRANSModelFactory> m_ebPatchRANSModelFactory;

  bool m_isDxSet;
  bool m_isDBLSet;
  bool m_isEBLGSet;
  bool m_isQuadCFISet;
  bool m_isEBISLSet;
  RealVect m_dx;
  int m_nGhost;
  int m_nComp;
  int m_nPrim;
  int m_nEqn;
  Vector<string> m_stateNames;
  Vector<string> m_primNames;
  bool m_hasCoarser, m_hasFiner;

  const EBLevelGrid* m_eblgPtr;

  LevelData<EBCellFAB>         m_stateOld;
  LevelData<EBCellFAB>         m_stateNew;

  LevelData<EBCellFAB>         m_redisRHS;
  LevelData<BaseIVFAB<Real> >  m_massDiff;
  LayoutData<IntVectSet>       m_sets;

  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  const RefCountedPtr<EBQuadCFInterp>*                          m_quadCFIPtr;

  EBLevelRedist                                                 m_ebLevelRedist;
  EBCoarseAverage                                               m_ebCoarseAverage;
  EBFluxRegister                                                m_divFFluxRegister;
  EBFluxRegister                                                m_scalarFluxRegister;
  EBPWLFineInterp                                               m_ebFineInterp;
  EBLevelRANS                                                   m_ebLevelRANS;

   //! Initial values of conserved quantities.
  Real m_originalMass;

  const LevelData<EBFluxFAB>*                      m_advVelPtr;
  const LayoutData< Vector< BaseIVFAB<Real> * > >* m_coveredAdvVelLoPtr;
  const LayoutData< Vector< BaseIVFAB<Real> * > >* m_coveredAdvVelHiPtr;

  LevelData<EBCellFAB> m_normalVelOld;
  LevelData<EBCellFAB> m_normalVelNew;
  bool m_oldNormalVelSet;

  // diffusion stuff
  // coefficients: aco*dS/dt = source + del.(bco*delS)
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        m_bco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bcoIrreg;

  Vector<RefCountedPtr<BaseDomainBCFactory> >         m_domBC;
  Vector<RefCountedPtr<BaseEBBCFactory> >             m_ebBC;

  // solvers:
  // AMRMG stuff:
  static Vector<RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > >  s_diffuseOpFactory;
  static Vector<RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > > > s_diffuseSolver;

  // bottom solver
//  static BiCGStabSolver<LevelData<EBCellFAB> > s_botSolver;
  static EBSimpleSolver                          s_botSolver;

  // temp integrators
  static Vector<RefCountedPtr<EBLevelBackwardEuler> >  s_diffuseIntegratorBE;
  static Vector<RefCountedPtr<EBLevelTGA> >            s_diffuseIntegratorTGA;

private():

  EBAMRRANS* getCoarserLevel() const;
  EBAMRRANS* getFinerLevel() const;

  void levelSetup();

  //weak construction is an
  //evil that must be rooted out.
  //verily, it is a character flaw.
  EBAMRRANS()
  {
    MayDay::Error("invalid operator");
  }
  //disallowed for all the usual reasons
  void operator=(const EBAMRRANS& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRRANS(const EBAMRRANS& a_input)
  {
    MayDay::Error("invalid operator");
  }
  
};

#include "NamespaceFooter.H"
#endif
