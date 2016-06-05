#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRTRANSPORT_H_
#define _EBAMRTRANSPORT_H_

#include "EBAMRLevel.H"
#include "EBLevelTransport.H"
#include "EBPatchTransport.H"
#include "EBAMRTransportParams.H"

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

// transport equation solver
// driven by external AMR or EBAMR class
class EBAMRTransport: public EBAMRLevel 
{
public:

  virtual ~EBAMRTransport();

// a_externalDriver = true is the grids, dt, advection vel, etc are set by external driver
  EBAMRTransport(const EBAMRTransportParams& a_params,
                 const RefCountedPtr<EBPatchGodunovFactory>& a_godFactory,
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

  virtual void initialGrid(const Vector<Box>& a_new_grids);

  virtual void initialData();

  virtual void postInitialize();

  virtual Real computeDt();

  virtual Real computeInitialDt();

  virtual void assignDBLPtr(const DisjointBoxLayout* a_dblPtr);

  virtual void assignEBLGPtr(const EBLevelGrid* a_eblgPtr);

  virtual void assignEBISLPtr(const EBISLayout* a_ebislPtr);

  virtual void assignQuadCFIPtr(const RefCountedPtr<EBQuadCFInterp>* a_quadCFIPtr);

  virtual void syncWithFineLevel();

  //! Sums the conserved variable identified by \a_ivar over the entire grid level, placing 
  //! the result in \a a_sumcons.
  //! \param a_sumcons The value of the conserved variable summed over the grid level.
  //! \param a_ivar The index of the conserved variable within its FAB.
  void sumConserved(Real& a_sumcons,
                    const int& a_ivar) const; 

  virtual void assignAdvVelPtr(const LevelData<EBFluxFAB>* a_advVelPtr,
                               const LayoutData< Vector< BaseIVFAB<Real> * > >* a_coveredAdvVelLoPtr,
                               const LayoutData< Vector< BaseIVFAB<Real> * > >* a_coveredAdvVelHiPtr);

  void defineSolvers();

  void postInitialGrid(const bool a_restart);

  virtual LevelData<EBCellFAB>& getConsStateNew();

  virtual Vector<string> getConsNames();

  virtual int numConserved();  

  virtual void setNormalVelOld(const LevelData<EBCellFAB>& a_normalVel);

  virtual void setNormalVelOld(const LevelData<EBFluxFAB>& a_advVel,
                               const LayoutData< Vector< BaseIVFAB<Real> * > >& a_coveredAdvVelLo,
                               const LayoutData< Vector< BaseIVFAB<Real> * > >& a_coveredAdvVelHi);

  #ifdef CH_USE_HDF5
  // ---------------------------------------------------
  // HDF5 routines for plotting and writing checkpoints.
  // ---------------------------------------------------

//  static int s_NewPlotFile;
  static bool s_solversDefined;
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
 
  void defineFactories(bool a_atHalfTime);

  void postTimeStepRefluxRedistDance();

  void implicitReflux();

  void refluxRedistInteraction();

  void resetWeights();

  // This performs the redistribution of a quantity (mass, momentum, or energy) in m_stateNew. 
  // The quantity is determined by a_interv. If a_putIntoState is set, the redistribution is 
  // done EXPLICITLY, and the redistributed quantity is placed directly into m_stateNew.
  // Otherwise, the quantity is accumulated into m_redisRHS for the later implicit 
  // redistribution.
  void refluxRHS();
  void explicitReflux(const Interval& a_interv);
  void coarseFineRedistribution(const Interval& a_interv);

  void getHalfState(LevelData<EBCellFAB>& a_stateInt);

  // inviscid advance
  void explicitAdvance(const LevelData<EBCellFAB>& a_divergeF);

  //(alpha I - dt Ls) delta = dt*Dr(Fs)
  void getRefluxDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                       Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeS,
                       int baseLev, int finestLev, Real baseDt); 

  //scalar += dt*rho*Cv*(deltaT)
  void incrScalarByDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                          int baseLev, int finestLev);
  void incrScalarByDeltaS(LevelData<EBCellFAB>& a_state,
                          LevelData<EBCellFAB>& a_deltaScalar);  

  void fluxDivergence(  LevelData<EBCellFAB>& a_divergeF);

  // advances the solution using the hyperbolic terms.
  //(\kappa \rho^\npo I - \frac{\kappa \dt}{2} L^m)\ubar^\npo = \kappa \rho^\npo \ubold^*
  void  getUStar(LevelData<EBCellFAB>      & a_UStar,
                 const LevelData<EBCellFAB>& a_UN,
                 const LevelData<EBCellFAB>& a_divergef);

  //if doNormalization = true, averages with neighboring cells
  //otherwise returns kappa * source
  void explicitHyperbolicSource(LevelData<EBCellFAB>&       a_diffSource,
                                const LevelData<EBCellFAB>& a_state,
                                bool a_doNormalization );

  //this one is explicit and also does reweighting of the 
  //redistribution object in the case of mass weighting
  void hyperbolicRedistribution(LevelData<EBCellFAB>& a_state); 

  void getDivDGradS(LevelData<EBCellFAB>& a_divDGradS,
                    LevelData<EBCellFAB>& a_Ustar);

  /// set  output to volfrac*((del dot (kappa grad S))
  // if there is viscous dissipation, need to add here
  void kappaDiffusiveSource(LevelData<EBCellFAB>& a_kappaDiffSource,
                            const LevelData<EBCellFAB>& a_state,
                            const LevelData<EBCellFAB>* a_stateCoar);

  //this includes news so need to call delete
  void
  getCoarserState(LevelData<EBCellFAB>* & a_stateCoar); 

  void coarseFineIncrement();

  void hyperbolicSource(LevelData<EBCellFAB>& a_source);

  //put ustar into state
  void finalAdvance(LevelData<EBCellFAB>& a_Ustar); 

  void setBCs();

  bool m_externalDriver;
  EBAMRTransportParams m_params;
  RefCountedPtr<EBPatchGodunovFactory> m_ebPatchGodunovFactory;

  bool m_isDxSet;
  bool m_isDBLSet;
  bool m_isEBLGSet;
  bool m_isQuadCFISet;
  bool m_isEBISLSet;
  RealVect m_dx;
  int m_nGhost;
  int m_nComp;
  int m_nPrim;
  Vector<string> m_stateNames;
  Vector<string> m_primNames;
  bool m_hasCoarser, m_hasFiner;

//  const DisjointBoxLayout* m_dblPtr;
  const EBLevelGrid* m_eblgPtr;
//  const RefCountedPtr<EBLevelGrid> m_eblgPtr;
//  EBLevelGrid        m_eblgOld;
//  const EBISLayout* m_ebislPtr;
  RefCountedPtr<EBPatchTransport> m_ebPatchGodunov;

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
  EBLevelTransport                                              m_ebLevelGodunov;

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
  RefCountedPtr<LevelData<EBCellFAB> >        m_aco;
  RefCountedPtr<LevelData<EBFluxFAB> >        m_bco;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_bcoIrreg;

  RefCountedPtr<BaseDomainBCFactory>          m_domBC;
  RefCountedPtr<BaseEBBCFactory>              m_ebBC;

  // solvers:
  // AMRMG stuff:
  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >  s_diffuseOpFactory;
  static RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >  s_diffuseSolver;
  
  // bottom solver
//  static BiCGStabSolver<LevelData<EBCellFAB> > s_botSolver;
  static EBSimpleSolver                        s_botSolver;  

  // temp integrators
  static RefCountedPtr<EBLevelBackwardEuler>   s_diffuseIntegratorBE;
  static RefCountedPtr<EBLevelTGA>             s_diffuseIntegratorTGA;  

  void clearSolvers();  

  void fillCoefficients(const LevelData<EBCellFAB>& a_state);

  void setDiffusionCoefficients(LevelData<EBCellFAB>& a_stateCell, 
                                LevelData<EBFluxFAB>& a_stateFace);

  int getFinestLevel();
 
private:
  
  EBAMRTransport* getCoarserLevel() const;  
  EBAMRTransport* getFinerLevel() const;

  void levelSetup();
  
  //weak construction is an
  //evil that must be rooted out.
  //verily, it is a character flaw.
  EBAMRTransport()
  {
    MayDay::Error("invalid operator");
  }
  //disallowed for all the usual reasons
  void operator=(const EBAMRTransport& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRTransport(const EBAMRTransport& a_input)
  {
    MayDay::Error("invalid operator");
  }  

};

#include "NamespaceFooter.H"
#endif
