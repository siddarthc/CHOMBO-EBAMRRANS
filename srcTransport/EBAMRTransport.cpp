#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "EBAMRTransport.H"
#include "parstream.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "DebugOut.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "SPMD.H"
#include "EBLoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "AMRMultiGrid.H"
#include "EBAMRIO.H"
#include "BaseIVFactory.H"
#include "EBViscousTensorOpFactory.H"
#include "EBConductivityOpFactory.H"
#include "EBAMRPoissonOpFactory.H"
#include "KappaSquareNormal.H"

#include "EBAMRLevel.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
//#include "EBLGIntegrator.H"
#include "EBLevelDataOps.H"
#include "EBBackwardEuler.H"
#include "AMR.H"
#include "NeumannConductivityDomainBC.H"
#include "NeumannConductivityEBBC.H"

#include "NamespaceHeader.H"

RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > EBAMRTransport::s_diffuseOpFactory        =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();

RefCountedPtr< EBLevelBackwardEuler>                     EBAMRTransport::s_diffuseIntegratorBE    =  RefCountedPtr< EBLevelBackwardEuler>();

RefCountedPtr< EBLevelTGA>                               EBAMRTransport::s_diffuseIntegratorTGA   =  RefCountedPtr< EBLevelTGA>();

RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > > EBAMRTransport::s_diffuseSolver          =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();

//BiCGStabSolver<LevelData<EBCellFAB> >                    EBAMRTransport::s_botSolver          = BiCGStabSolver<LevelData<EBCellFAB> >();

EBSimpleSolver                 EBAMRTransport::s_botSolver;
bool EBAMRTransport::s_solversDefined = false;
/******************/
EBAMRTransport::
~EBAMRTransport()
{
//  delete m_dblPtr;
//  if (m_eblgPtr != NULL) delete m_eblgPtr;

//  delete m_ebislPtr;
//  if (m_quadCFIPtr != NULL) delete m_quadCFIPtr;

// these deletes are called in driver class
//  if (m_advVelPtr !=NULL) delete m_advVelPtr;

//  if (m_coveredAdvVelLoPtr != NULL) delete m_coveredAdvVelLoPtr;

//  if (m_coveredAdvVelHiPtr != NULL) delete m_coveredAdvVelHiPtr;

// deletes in base class
//  if (m_dtPtr != NULL) delete m_dtPtr;

//  if (m_timePtr != NULL) delete m_timePtr;

  if ((m_level == 0) && m_params.m_doDiffusion)
    {
      clearSolvers();
    }
}
/****************/
void EBAMRTransport::clearSolvers()
{
  s_diffuseOpFactory = RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
  s_diffuseIntegratorBE = RefCountedPtr< EBLevelBackwardEuler>();
  s_diffuseIntegratorTGA = RefCountedPtr< EBLevelTGA>();
  s_diffuseSolver = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >();
}
/******************/
/******************/
EBAMRTransport::
EBAMRTransport(const EBAMRTransportParams& a_params,
               const RefCountedPtr<EBPatchGodunovFactory>& a_godFactory,
               bool  a_externalDriver):
    m_params(a_params),
    m_ebPatchGodunovFactory(a_godFactory),
    m_externalDriver(a_externalDriver)
{
  m_isDxSet = false;
  m_isDBLSet = false;
  m_isEBLGSet = false;
  m_isQuadCFISet = false;
  m_isEBISLSet = false;
//  m_dblPtr = NULL;
  m_eblgPtr = NULL;
//  m_ebislPtr = NULL;

  m_oldNormalVelSet = false;

}
/******************/ 
/******************/
void EBAMRTransport::
assignDx(RealVect a_dx)
{
  m_isDxSet = true;
  m_dx = a_dx;
}
/******************/
/******************/
void EBAMRTransport::
assignDBLPtr(const DisjointBoxLayout* a_dblPtr)
{
//  m_isDBLSet = true;
//  m_dblPtr = a_dblPtr;
}
/******************/
/******************/
void EBAMRTransport::
assignEBLGPtr(const EBLevelGrid* a_eblgPtr)
{
  m_isEBLGSet = true;
  m_eblgPtr = a_eblgPtr;
}
/******************/
/******************/
void EBAMRTransport::
assignQuadCFIPtr(const RefCountedPtr<EBQuadCFInterp>* a_quadCFIPtr)
{
  m_isQuadCFISet = true;
  m_quadCFIPtr = a_quadCFIPtr;
}
/******************/
/******************/
void EBAMRTransport::
assignEBISLPtr(const EBISLayout* a_ebislPtr)
{
//  m_isEBISLSet = true;
//  m_ebislPtr = a_ebislPtr;
}
/******************/
/******************/
void EBAMRTransport::
define(EBAMRLevel*            a_coarser_level_ptr,
       const ProblemDomain& a_problem_domain,
       int                  a_level,
       int                  a_ref_ratio)
{
  if (m_params.m_verbosity >=3)
    {
      pout() << "EBAMRTransport::define, level=" << a_level << endl;
    }

  EBPatchGodunov::useConservativeSource(true);
  m_isDefined = true;
  EBAMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);

  if (a_coarser_level_ptr != NULL)
    { 
      EBAMRTransport* amrg_ptr = 
        dynamic_cast<EBAMRTransport*>(a_coarser_level_ptr);
      if (amrg_ptr == NULL)
        {
          pout() << "EBAMRG::define:cannot cast  to derived class"
                 << endl;
          MayDay::Error();
        }
      m_params = amrg_ptr->m_params;
    }

  if (m_externalDriver)
    {
      CH_assert(m_isDxSet );
    }
  else 
    {
      pout() << "EBAMRTransport class is not setup to drive itself" << endl;
      MayDay::Error();
    }

  m_nGhost = 4; // hmm...
  
  m_ebPatchGodunov = RefCountedPtr<EBPatchTransport>();
  m_ebPatchGodunov = RefCountedPtr<EBPatchGodunov>(m_ebPatchGodunovFactory->create());
  m_ebPatchGodunov->define(m_problem_domain, m_dx);

  m_nComp      = m_ebPatchGodunov->numConserved();
  m_nPrim      = m_ebPatchGodunov->numPrimitives();
  m_stateNames = m_ebPatchGodunov->stateNames();
  m_primNames  = m_ebPatchGodunov->primNames();
  m_ref_ratio  = a_ref_ratio;
}
/******************/
int EBAMRTransport::numConserved()
{
  return m_ebPatchGodunov->numConserved();
}
/******************/
LevelData<EBCellFAB>& EBAMRTransport::getConsStateNew()
{
  return m_stateNew;
}
/******************/ 
Vector<string> EBAMRTransport::getConsNames()
{
  return m_stateNames;
}
/******************/
void EBAMRTransport::getHalfState(LevelData<EBCellFAB>& a_stateInt)
{
  //interpolate state to n+1/2
  Real told = 0; Real tnew = 1; Real time = 0.5;
  EBArith::timeInterpolate(a_stateInt, m_stateOld, m_stateNew,
                           m_eblgPtr->getDBL(), time, told, tnew);
}
/******************/
Real EBAMRTransport::advance()
{
  pout() << "advancing EBAMRTransport with dt = " << m_dt << endl;
  
  EBPatchGodunov::s_whichLev = m_level;

  if(m_params.m_variableCoeff && (m_level== 0) && m_params.m_doDiffusion) defineSolvers(); 

  m_stateNew.copyTo(m_stateNew.interval(), m_stateOld, m_stateOld.interval());

  EBCellFactory fact(m_eblgPtr->getEBISL());
  LevelData<EBCellFAB>      divergeF(m_eblgPtr->getDBL(),m_nComp          , 4*IntVect::Unit, fact);

  pout() << "EBAMRTransport::getting diverence of flux" << endl;

  fluxDivergence( divergeF);

  if(!m_params.m_doDiffusion)
    {
      explicitAdvance(divergeF);
    }

  else
    {
//      pout() << "EBAMRTransport::diffusion is not set" << endl;
      // diffusion fance
      pout() << "EBAMRTransport::step1: getting Ustar (U^n + explicit contribution)" << endl;

      //U* = Un - dt*L^c(U^n)
      LevelData<EBCellFAB>       UStar(m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);
      getUStar(UStar, m_stateOld, divergeF);

      pout() << "step 1.5: doing explicit redistribution into ustar of hyperbolic mass difference" <<endl;
      //redistribute stuff in hyperbolic registers.
      hyperbolicRedistribution(UStar);

      ////debug
      //LevelData<EBCellFAB> consAndPrimStar;
      //fillConsAndPrim(consAndPrimStar, UStar, 1);
      ////end debug

      // for viscous dissipation, need to do the divSigma stuff.. neglecting for now
      if (m_params.m_doViscDissipation)
        {
          MayDay::Error("EBAMRTransport::viscous dissipation is not set");
        }     

      // now add the div(kappaGradS)
      pout() << "step2: solve for diffusion term and add to transport equation" <<endl;
      LevelData<EBCellFAB> dtDivDGradS(m_eblgPtr->getDBL(),       1, 4*IntVect::Unit, fact);
      getDivDGradS(dtDivDGradS, UStar);

      pout() << "putting state into m_statenew and flooring" << endl;
      finalAdvance(UStar);
    }
  
  if(m_params.m_checkMaxMin)
    {
      CH_TIME("max_min_check");
      pout() << "EBAMRTransport::after advance max mins for data for level " << m_level << endl;
      for(int icomp = 0; icomp < m_stateNew.nComp(); icomp++)
        {
          Real maxVal, minVal;
          EBLevelDataOps::getMaxMin(maxVal, minVal, m_stateNew, icomp);
          pout() << " "
                 << setprecision(4)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 <<  "comp = " <<  icomp
                 << ", max = " << maxVal
                 << ", min = " << minVal << endl;
        }
    } 

  return m_dt;
}
/******************/
/******************/
//U* = Un + dt*L^c(U^n) + (dt/2)*L^v(U^n)
void EBAMRTransport::
getUStar(LevelData<EBCellFAB>      & a_UStar,
         const LevelData<EBCellFAB>& a_UN,
         const LevelData<EBCellFAB>& a_divergeF)
{
  EBCellFactory fact(m_eblgPtr->getEBISL());
  //advance everything explicitly
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB  dtLcU(m_eblgPtr->getEBISL()[dit()], m_eblgPtr->getDBL().get(dit()), m_nComp);
      dtLcU.setVal(0.);
      dtLcU += a_divergeF[dit()];
      dtLcU *= m_dt;

      a_UStar[dit()].setVal(0.);
      a_UStar[dit()] += a_UN[dit()];
      a_UStar[dit()] -= dtLcU;
    }

  m_ebLevelGodunov.floorConserved(a_UStar, m_time, m_dt);
}
/******************/
/******************/
void EBAMRTransport::
fluxDivergence(  LevelData<EBCellFAB>& a_divergeF)
{
  CH_assert(a_divergeF.nComp() == m_nComp);
  CH_assert(m_externalDriver);
  CH_assert(m_oldNormalVelSet);
  CH_assert(m_oldTimeSet);
  CH_assert(m_advVelPtr != NULL);
  CH_assert(m_coveredAdvVelLoPtr != NULL);
  CH_assert(m_coveredAdvVelHiPtr != NULL);

  pout() << "EBAMRTransport::taking flux divergence on level " << m_level << " with dt = " << m_dt << endl;
  IntVect ivGhost = 4*IntVect::Unit;
  EBCellFactory fact(m_eblgPtr->getEBISL());
  LevelData<EBCellFAB> source(m_eblgPtr->getDBL(),m_nComp, 4*IntVect::Unit, fact);
  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;
  const LevelData<EBCellFAB> ldv;
//  const LevelData<EBFluxFAB> ldf;
  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;
  const LevelData<EBCellFAB>* coarNormVelOld = &ldv;
  const LevelData<EBCellFAB>* coarNormVelNew = &ldv;

  Real told = 0.0;
  Real tnew = 0.0;

  m_ebLevelGodunov.computeNormalVel(m_normalVelNew, *m_advVelPtr, *m_coveredAdvVelLoPtr, *m_coveredAdvVelHiPtr);

  if(m_hasCoarser)
    {
      EBAMRTransport* coarPtr = getCoarserLevel();
      //recall that my flux register goes between my
      //level and the next finer level
      coarFR = &coarPtr->m_divFFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarNormVelOld = &coarPtr->m_normalVelOld;    

      // just to be consistent with the subcycling world
 
      coarDataNew = &coarPtr->m_stateNew;
      coarNormVelNew = &coarPtr->m_normalVelNew;
      tnew = coarPtr->m_time;
      told = tnew - coarPtr->m_dt;

      // At this level, m_time theoretically can not be smaller than its coarser level time (told); 
      // if it happends due to the precision/floating error, then we resign m_time to be told. 
      // One alternative way is to ensure a tolerance check. However, this current fix is intend to 
      // avoid further error accumulation (based on the choice of tolerance). 
      // Ideally, one would want to have a "generic" fix in either EBAMRNoSubcycle.cpp or postTimeStep().
      if (m_time < told)
        {
          m_time = told;
        }
    }
  if(m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_divFFluxRegister;
    }

  hyperbolicSource(source);
  // flux register manipulation happens in levelgodunov

  m_ebLevelGodunov.divergeF(a_divergeF,
                            m_massDiff,
                            *fineFR,
                            *coarFR,
                            m_stateOld,
                            m_normalVelNew,
                            *m_advVelPtr,
                            *m_coveredAdvVelLoPtr,
                            *m_coveredAdvVelHiPtr,
                            source,
                            *coarDataOld,
                            *coarDataNew,
                            *coarNormVelOld,
                            *coarNormVelNew,
                            m_time,
                            told,
                            tnew,
                            m_dt);

  coarseFineIncrement();
}
/******************/
/******************/
void EBAMRTransport::
hyperbolicSource(LevelData<EBCellFAB>&       a_source)
{
  bool zeroHyperbolicSource = false;
  
  // begin debug 
//  zeroHyperbolicSource = true;
  //set source to zero
  if(zeroHyperbolicSource)
    {
      pout() << "zeroing out hyperbolic source"  << endl;
      EBLevelDataOps::setToZero(a_source);
    }
  //end debug
  else
    {
      EBPatchGodunov::useConservativeSource(true);
      //this is important because sometimes there is no diffusion
      EBLevelDataOps::setVal(a_source, 0.0);
      if(m_params.m_doDiffusion)
        {
//          pout() <<"EBAMRTransport::diffusion not set "<< endl;
          EBCellFactory fact(m_eblgPtr->getEBISL());
          LevelData<EBCellFAB>  diffSrc(m_eblgPtr->getDBL(), m_nComp,4*IntVect::Unit, fact);
          explicitHyperbolicSource(diffSrc, m_stateOld, true);
          
          Interval interv(0, 0);
          for(DataIterator dit = a_source.dataIterator(); dit.ok(); ++dit)
            {
              const Box& region = m_eblgPtr->getDBL().get(dit());
              a_source[dit()].copy(region, interv, region, diffSrc[dit()], interv);

            }
          a_source.exchange();
        } 
     }
}
/******************/
/******************/
void EBAMRTransport::
explicitHyperbolicSource(LevelData<EBCellFAB>&       a_diffSource,
                         const LevelData<EBCellFAB>& a_state,
                         bool a_doNormalization)
{
  pout() << "computing explicit hyperbolic source" << endl;
  EBCellFactory fact(m_eblgPtr->getEBISL());
//  LevelData<EBCellFAB>  state(m_eblgPtr->getDBL(), m_nComp,4*IntVect::Unit, fact);  
  
  LevelData<EBCellFAB>* stateCoar = NULL;
  getCoarserState(stateCoar);

  kappaDiffusiveSource(a_diffSource, a_state, stateCoar);

  if (m_hasCoarser)
    {
      delete stateCoar;
    } 
  
  if (a_doNormalization)
    {
      //finally all these things are multiplied by kappa so we need to make 
      //them normalized 
      KappaSquareNormal normalizinator(*m_eblgPtr);
      normalizinator(a_diffSource);
    }  
}
//---------------------------------------------------------------------------------------
/// set  output to (del dot (kappa grad S))  
// if there is viscous dissipation, need to add here
void
EBAMRTransport::
kappaDiffusiveSource(LevelData<EBCellFAB>& a_kappaDiffSource,
                     const LevelData<EBCellFAB>& a_state,
                     const LevelData<EBCellFAB>* a_stateCoar)
{
  //set the a coefficients to time n
  Interval srcInt(0, 0);
  Interval dstInt(0, 0);

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblgPtr->getDBL().get(dit());
      (*m_aco)[dit()].setVal(m_params.m_alphaCoeff);
    }

 // Compute the diffusion term = volfrac(div kappa grad S)
  Real alpha = 0;   Real beta = 1; //want just the div(flux) part of the operator
  s_diffuseIntegratorBE->applyOperator(a_kappaDiffSource, a_state,  a_stateCoar, m_level, alpha, beta, true); 

  if (m_params.m_doViscDissipation)
    {
      MayDay::Error("EBAMRTransport:: viscous dissipation is not set");
    }

  /**/
}
/******************/
/******************/
// this includes news so need to call delete
void EBAMRTransport::
getCoarserState(LevelData<EBCellFAB>* & a_stateCoar)
{
  LevelData<EBCellFAB>* stateCoar = NULL;
  if(m_hasCoarser)
    {
      EBAMRTransport* coarPtr = getCoarserLevel();
      EBCellFactory factCoar(coarPtr->m_eblgPtr->getEBISL());
      stateCoar = new LevelData<EBCellFAB>(coarPtr->m_eblgPtr->getDBL(), m_nComp,4*IntVect::Unit, factCoar);
      LevelData<EBCellFAB>    stateCoarOld(coarPtr->m_eblgPtr->getDBL(), m_nComp,4*IntVect::Unit, factCoar);
      LevelData<EBCellFAB>    stateCoarNew(coarPtr->m_eblgPtr->getDBL(), m_nComp,4*IntVect::Unit, factCoar);

      coarPtr->m_stateOld.copyTo(stateCoarOld);
      coarPtr->m_stateNew.copyTo(stateCoarNew);

      Real tCoarNew = coarPtr->m_time;
      Real tCoarOld = tCoarNew - coarPtr->m_dt;

      //interpolate coarse solution to fine time
      EBArith::timeInterpolate(*stateCoar, stateCoarOld, stateCoarNew, coarPtr->m_eblgPtr->getDBL(), m_time, tCoarOld, tCoarNew);
    }

  a_stateCoar = stateCoar;
}
/******************/
/******************/
void EBAMRTransport::
coarseFineIncrement()
{
  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level
  {
    CH_TIME("coarse-fine guano");
    if(m_hasCoarser && (!s_noEBCF))
      {
        for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
          {
            m_ebFineToCoarRedist.increment(m_massDiff[dit()], dit(), interv);
          }
      }

    //initialize redistribution register between this level
    // and next finer level
    //this includes re-redistribution registers
    if(m_hasFiner  && (!s_noEBCF))
      {
        m_ebCoarToFineRedist.setToZero();
        m_ebCoarToCoarRedist.setToZero();
        for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
          {
            BaseIVFAB<Real>& massDiffFAB = m_massDiff[dit()];
            m_ebCoarToFineRedist.increment(massDiffFAB, dit(), interv);
            m_ebCoarToCoarRedist.increment(massDiffFAB, dit(), interv);
          }
      }
  }
}
/******************/
/******************/
int EBAMRTransport::
getFinestLevel()
{
  int imax = 0;
  Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
  for(int ilev = 0; ilev < hierarchy.size(); ilev++)
    {
      EBAMRTransport* transportLevel = dynamic_cast<EBAMRTransport*>(hierarchy[ilev]);
      int numBoxes = 0;
      if(transportLevel->m_eblgPtr->getDBL().isClosed()) numBoxes = transportLevel->m_eblgPtr->getDBL().size();
      if(numBoxes > 0)
        imax++;
      else
        break;
    }
  return imax-1;
}
/******************/
/******************/
void EBAMRTransport::
explicitAdvance(const LevelData<EBCellFAB>& a_divergeF)
{
  CH_assert(!m_params.m_doDiffusion);
  //advance everything explicitly
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB dtDivergeF(m_eblgPtr->getEBISL()[dit()], m_eblgPtr->getDBL().get(dit()), a_divergeF.nComp());
      dtDivergeF.setVal(0.);
      dtDivergeF += a_divergeF[dit()];
      dtDivergeF *= m_dt;
      m_stateNew[dit()] -= dtDivergeF;
    }
  hyperbolicRedistribution(m_stateNew);

  m_ebLevelGodunov.floorPrimitives(m_stateNew, m_timeOld, m_dtOld);
}
/******************/
/******************/
void EBAMRTransport::
hyperbolicRedistribution(LevelData<EBCellFAB>& a_state)
{
  //explicit redistribution here because these are the hyperbolic terms
  m_ebLevelRedist.setToZero();
  if(m_params.m_useMassRedist)
    {
      //if use mass weighting, need to
      //fix weights of redistribution object
      int densityIndex = m_ebPatchGodunov->densityIndex();
      a_state.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(a_state, densityIndex);
    }

  Interval consInterv(0, m_nComp-1);
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), consInterv);
    }
  //
  m_ebLevelRedist.redistribute(a_state, consInterv);

  m_ebLevelRedist.setToZero();
}
/******************/
/******************/
void EBAMRTransport::
getDivDGradS(LevelData<EBCellFAB>& a_dtDivDGradS,
             LevelData<EBCellFAB>& a_UStar)
{
  EBCellFactory fact(m_eblgPtr->getEBISL());
  LevelData<EBCellFAB>  Sold(m_eblgPtr->getDBL(),        1,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>  Snew(m_eblgPtr->getDBL(),        1,4*IntVect::Unit, fact);

  a_UStar.copyTo(Sold);

  EBFluxRegister*       coarFRPtr = NULL;
  EBFluxRegister*       fineFRPtr = NULL;
  LevelData<EBCellFAB>* SCoarOldPtr = NULL;
  LevelData<EBCellFAB>* SCoarNewPtr = NULL;
  Real tCoarOld= 0;
  Real tCoarNew= 0;

  if(m_hasCoarser)
    {
      // Get coarser data if we have it.
      EBAMRTransport* coarTransport = getCoarserLevel();
      coarFRPtr = &coarTransport->m_scalarFluxRegister;
      tCoarNew = coarTransport->m_time;
      tCoarOld = tCoarNew - coarTransport->m_dt;

      const EBLevelGrid& ceblg = *coarTransport->m_eblgPtr;
      EBCellFactory cfact(ceblg.getEBISL());
      SCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, 4*IntVect::Unit, cfact);
      SCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, 4*IntVect::Unit, cfact);

      coarTransport->m_stateOld.copyTo(*SCoarOldPtr);
      coarTransport->m_stateNew.copyTo(*SCoarNewPtr);
    }

  if(m_hasFiner)
    {

      // Get finer flux registers if we have them.

      fineFRPtr = &m_scalarFluxRegister;
    }
  // Set the a coefficients for the thermal conduction equation to Cv * rho, 
  // specifying both the old and new values so that the density gets linearly
  // interpolated.
  Interval srcInt(0, 0);
  Interval dstInt(0, 0);

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblgPtr->getDBL().get(dit());
      (*m_aco)[dit()].setVal(m_params.m_alphaCoeff);
    }

  LevelData<EBCellFAB> rhs(m_eblgPtr->getDBL(), 1,  m_nGhost*IntVect::Unit, fact);
   EBLevelDataOps::setToZero(rhs);
  Interval srcComp(0,0);
  Interval dstComp(0,0);

  m_redisRHS.copyTo(srcComp, rhs, dstComp);
  EBLevelDataOps::setToZero(m_redisRHS);

  if(m_params.m_backwardEuler)
    {
      s_diffuseIntegratorBE->updateSoln(Snew, Sold, rhs,  fineFRPtr, coarFRPtr,
                                     SCoarOldPtr, SCoarNewPtr, m_time, tCoarOld,
                                     tCoarNew, m_dt/2.,m_level, true);
    }
  else
    {
      s_diffuseIntegratorTGA->updateSoln(Snew, Sold, rhs,  fineFRPtr, coarFRPtr,
                                     SCoarOldPtr, SCoarNewPtr, m_time, tCoarOld,
                                     tCoarNew, m_dt/2.,  m_level, true);

    }

  EBLevelDataOps::scale(Sold, -1.0);
  EBLevelDataOps::scale(Snew,  1.0);

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_dtDivDGradS[dit()].setVal(0.);
      //sets divDgradS = ((Snew-Sold)/dt)*dt
      a_dtDivDGradS[dit()] += Sold[dit()]; //neg sign is in the scale call above
      a_dtDivDGradS[dit()] += Snew[dit()];

      //this makes divDgradS = (alpha (Snew-Sold)/dt)*dt
      a_dtDivDGradS[dit()] *= (*m_aco)[dit()];

      int isrc = 0; int idst = 0; int inco = 1;
      a_UStar[dit()].plus(a_dtDivDGradS[dit()], isrc, idst, inco);
    }

  if(m_hasCoarser)
    {
      delete SCoarOldPtr;
      delete SCoarNewPtr;
    }
}
/******************/
/******************/
void EBAMRTransport::
finalAdvance(LevelData<EBCellFAB>& a_Ustar)
{
  //set stateNew to udbst
  Interval comps(0, m_nComp-1);
  a_Ustar.copyTo(comps, m_stateNew, comps);

  m_ebLevelGodunov.floorConserved(m_stateNew, m_time, m_dt);   
} 
/******************/
/******************/
void EBAMRTransport::postTimeStep()
{
//  pout() << " in EBAMRTransport postTimeStep for level " << m_level << endl;

 // Average the finer grid levels to be consistent with this one.

  m_oldNormalVelSet = false;

  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRTransport* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }

 // this does the refluxing and redistribution evil dance
  postTimeStepRefluxRedistDance();
  m_ebLevelGodunov.floorPrimitives(m_stateNew, m_time, m_dt); 

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
      m_redisRHS[dit()].setVal(0.0);
    }
}
/******************/
/******************/
void
EBAMRTransport::postTimeStepRefluxRedistDance()
{
  resetWeights();
  refluxRedistInteraction();
  coarseFineRedistribution(Interval(0, m_nComp-1));
  explicitReflux(Interval(0, m_nComp-1));
  ParmParse pp;
  bool turn_off_reflux = false;
  pp.query("turn_off_implicit_reflux", turn_off_reflux);
  if(turn_off_reflux)
    {
      pout() << "implicit reflux turned off" << endl;
    }
  if(m_params.m_doDiffusion && !turn_off_reflux)
    {
      refluxRHS();
      //      defineSolvers();
      implicitReflux();
    }
}
/******************/
/******************/
 void EBAMRTransport::implicitReflux()
{
//  MayDay::Error("EBAMRTransport::implicitReflux is not set");
  pout() << "EBAMRTransport::implicit redist/reflux for level" << m_level << endl;
  Real crseTime = -1.0;
  Real timeEps = 1.0e-2*m_dt;
  if (m_level > 0) crseTime = m_coarser_level_ptr->time();

  if (m_level == 0 || (abs(crseTime - m_time) < timeEps))
    {
       int finestLev = getFinestLevel();
      //if I have not caught up to the next coarser level,
      //then I need to do implicit{reflux,redist} for my level and all
      //finer levels
      Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
      Vector<LevelData<EBCellFAB>* >    deltaScalar(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* > dtRefluxDivergeS(finestLev+1, NULL);
      for(int ilev = 0; ilev <= finestLev; ilev++)
        {
          EBAMRTransport* transportLevel = dynamic_cast<EBAMRTransport*>(hierarchy[ilev]);

          EBCellFactory fact(transportLevel->m_eblgPtr->getEBISL());
          deltaScalar[ilev]    = new LevelData<EBCellFAB>(transportLevel->m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);
          dtRefluxDivergeS[ilev] = new LevelData<EBCellFAB>(transportLevel->m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);

          //copy redistRHS to reflux divergence holders
          Interval srcComp, dstComp;
          srcComp = Interval(0, 0);
          dstComp = Interval(0, 0);
          transportLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeS[ilev], dstComp);

          EBLevelDataOps::scale(*dtRefluxDivergeS[ilev], transportLevel->m_dt);

          //these MATTER especially with subcycling
          EBLevelDataOps::setVal(*deltaScalar[ilev], 0.0);
        }

      pout() << "getting scalar reflux/redist increment"  << endl;
      int baseLev = Max(m_level-1, 0);
      Real baseDt = ((EBAMRTransport*)hierarchy[baseLev])->m_dt;

      //(alpha I - dt Ls) delta =dt* Dr(Fs)
      getRefluxDeltaS(deltaScalar, dtRefluxDivergeS, baseLev, finestLev, baseDt);

      //rho += dt* Ls(deltaS)
      incrScalarByDeltaS(deltaScalar, baseLev, finestLev);
    
      for(int ilev = 0; ilev <= finestLev; ilev++)
        {
          delete    deltaScalar[ilev];
          delete dtRefluxDivergeS[ilev];
        }

      //reset redistribution RHS to avoid double counting
      for(int ilev = m_level; ilev <= finestLev; ilev++)
        {
          EBAMRTransport* transportLevel = dynamic_cast<EBAMRTransport*>(hierarchy[ilev]);
          EBLevelDataOps::setToZero(transportLevel->m_redisRHS);
        }
    } 
}
/******************/
/******************/
//(alpha I - dt Ls) delta = dt*Dr(Fs)
void EBAMRTransport::
getRefluxDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeS,
                int a_baseLev, int a_finestLev, Real a_dtBase)
{
  Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for(int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRTransport* transportLevel = dynamic_cast<EBAMRTransport*>(hierarchy[ilev]);
      LevelData<EBCellFAB>& acoef = *transportLevel->m_aco;
      for(DataIterator dit = transportLevel->m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = transportLevel->m_eblgPtr->getDBL().get(dit());
          acoef[dit()].setVal(m_params.m_alphaCoeff);
        }
    }
  //set alpha to 1, beta = -dtbase
  s_diffuseIntegratorBE->resetSolverAlphaAndBeta(1.0, -a_dtBase);

  //solve equation (alpha I - dt Ls) delta = dt*Dr(Fs)
  //rhs already multiplied by dt
  //first true = zero phi
  //second true = force homogeneous bcs (solving for a delta S)
  s_diffuseSolver->solve(a_deltaScalar, a_dtRefluxDivergeS, a_finestLev, a_baseLev, true, true);
}
/******************/
/******************/
// scalar += dt*Ls(deltaS)
void EBAMRTransport::
incrScalarByDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                    int a_baseLev, int a_finestLev)
{
  Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRTransport* transportLevel   = dynamic_cast<EBAMRTransport*>(hierarchy[ilev]);
      transportLevel->incrScalarByDeltaS(transportLevel->m_stateNew, *a_deltaScalar[ilev]);
    }
} 
/******************/
/******************/
// scalar += dt*Ls(deltaS)
 void EBAMRTransport:: 
incrScalarByDeltaS(LevelData<EBCellFAB>& a_state,
                   LevelData<EBCellFAB>& a_deltaScalar)
{
  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      //add change in energy to the state
      {
        int isrc = 0; int idst = 0; int inco = 1;
        a_state[dit()].plus(a_deltaScalar[dit()], isrc, idst, inco);
      }
    } 
}
/******************/
/******************/
 void EBAMRTransport::refluxRedistInteraction()
{
  //the flux register must modify the redistribution
  //registers
  if(m_hasFiner  && (!s_noEBCF))
    {
      Real scale = -1.0/m_dx[0];
      Interval interv(0, m_nComp-1);
      m_divFFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist,
                                                 interv, scale);

      m_divFFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist,
                                                 interv, scale);
    } 
}
/******************/
/******************/
 void EBAMRTransport::resetWeights()
{
  //if use mass weighting, need to
  //fix weights of redistribution objects
  if(m_params.m_useMassRedist)
    {
      int densevar = m_ebPatchGodunov->densityIndex();
      m_stateNew.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(m_stateNew, densevar);

      if(m_hasCoarser && (!s_noEBCF))
        {
          EBAMRTransport* coarPtr = getCoarserLevel();
          coarPtr->m_stateNew.exchange(Interval(0, m_nComp-1));
          m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, densevar);
        }
      if(m_hasFiner && (!s_noEBCF))
        {
          m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
          m_ebCoarToCoarRedist.resetWeights(m_stateNew, densevar);
        }
    }
} 
/******************/
/******************/
void EBAMRTransport::explicitReflux(const Interval& a_interv)
{
  if (m_hasFiner)
    {
      Real scale = -1.0/m_dx[0];
      m_divFFluxRegister.reflux(m_stateNew, a_interv, scale);
      m_divFFluxRegister.setToZero();
    }
}
/******************/
/******************/
void EBAMRTransport::coarseFineRedistribution(const Interval& a_interv)
{
  if(m_hasCoarser)
    {
      if(m_params.m_doSmushing && (!s_noEBCF))
        {
          //redistibute to coarser level
          EBAMRTransport* coarPtr = getCoarserLevel();
          //put mass directly into state
          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, a_interv);
          m_ebFineToCoarRedist.setToZero();
        }
    }
  if (m_hasFiner)
    {
      EBAMRTransport* finePtr = getFinerLevel();
      if(m_params.m_doSmushing  && (!s_noEBCF))
        {
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, a_interv);
          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(         m_stateNew, a_interv);
          m_ebCoarToFineRedist.setToZero();
          m_ebCoarToCoarRedist.setToZero();
        }
    }
}
/******************/
/******************/
 void EBAMRTransport::refluxRHS()
{
//  MayDay::Error("EBAMRTransport::refluxRHS is not set");
  //this does the refluxing and redistribution evil dance
  Interval interv(0, 0);
  if (m_hasFiner)
    {
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];
      Interval interflux(0, 0);
      m_scalarFluxRegister.reflux(m_redisRHS, interv, interflux, scale);

      m_scalarFluxRegister.setToZero();
    } 
}
/******************/
/******************/
void EBAMRTransport::tagCells(IntVectSet& a_tags)
{
}
/******************/
/******************/
void EBAMRTransport::tagCellsInit(IntVectSet& a_tags)
{
}
/******************/
/******************/
void EBAMRTransport::regrid(const Vector<Box>& a_new_grids)
{
  // only for the case of being driven by external driver for now
  CH_assert(m_externalDriver);
  CH_assert(m_oldEBLGSet);
  CH_assert(m_newEBLGSet);
  CH_assert(m_isEBLGSet);  

  //first save old data
  // save data for later copy

  // eblgOld is eblg before regrid
  // it is assumed that regridding already happened in the driver class
  //      eblgNew is set through the pointer

  EBLevelGrid eblgOld = m_eblgOld;
  EBLevelGrid eblgNew = *m_eblgPtr;
  
  EBISLayout ebislOld = eblgOld.getEBISL();
  EBISLayout ebislNew = eblgNew.getEBISL();

  EBCellFactory factoryOld(ebislOld);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  LevelData<EBCellFAB> stateSaved(eblgOld.getDBL(), m_nComp, ivGhost, factoryOld);
  Interval interv(0,m_nComp-1);
  stateSaved.define(eblgOld.getDBL(), m_nComp, ivGhost, factoryOld);
  m_stateNew.copyTo(interv, stateSaved, interv);

  //four ghost cells.
  m_level_grids = a_new_grids;

  // set up data structures
  levelSetup();

  // interpolate to coarser level
  if (m_hasCoarser)
    {
      EBAMRTransport* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew,
                                 coarPtr->m_stateNew,
                                 interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
  m_stateNew.copyTo(interv,m_stateOld, interv);
}
/******************/
/******************/
void EBAMRTransport::
preRegrid(int                         a_base_level,
          const Vector<Vector<Box> >& a_new_grids)
{
   
}
/******************/
/******************/
void EBAMRTransport::postRegrid(int a_base_level)
{
  if(m_level == a_base_level) defineSolvers();
}
/******************/
/******************/
void EBAMRTransport::postInitialGrid(const bool a_restart)
{
  if(m_level == 0) defineSolvers();
}
/******************/
/******************/
void EBAMRTransport::defineSolvers()
{
  // need to set s_noEBCF to true if(tagAllIrreg). This is done in the driver class

  if (!m_externalDriver)
    {
      ParmParse pp;
      bool tagAllIrregular = false;
      if(pp.contains("tag_all_irregular"))
        {
          pp.get("tag_all_irregular", tagAllIrregular);
        }
      if(tagAllIrregular) s_noEBCF = true;
    }   

  EBConductivityOp::setForceNoEBCF( s_noEBCF);
  defineFactories(true);

  if (m_params.m_doDiffusion)
   {
//     MayDay::Error("diffusion set not set in EBAMRTransport");     
     Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);

      EBAMRTransport* coarsestLevel = dynamic_cast<EBAMRTransport*>(hierarchy[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblgPtr->getDomain();
      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRTransport* transportLevel = (EBAMRTransport*)(hierarchy[ilev]);
          eblgs       [ilev] = *transportLevel->m_eblgPtr;
          grids       [ilev] = transportLevel->m_eblgPtr->getDBL();
          refRat      [ilev] = transportLevel->m_ref_ratio;
        }

      //alpha, beta get replaced in tga solves
      IntVect  giv    = 4*IntVect::Unit;
      RealVect origin =  RealVect::Zero;
      int numSmooth, numMG, maxIter, mgverb;
      Real tolerance, hang, normThresh;
      ParmParse pp("scal_amrmultigrid");
      pp.get("num_smooth", numSmooth);
      pp.get("num_mg",     numMG);
      pp.get("hang_eps",   hang);
      pp.get("norm_thresh",normThresh);
      pp.get("tolerance",  tolerance);
      pp.get("max_iter",   maxIter);
      pp.get("verbosity",  mgverb);

      s_diffuseSolver = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());

      s_diffuseSolver->define(lev0Dom, *s_diffuseOpFactory, &s_botSolver, nlevels);

      s_diffuseSolver->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);

      s_diffuseSolver->m_verbosity = mgverb;
      Real bottomCushion = 1.0;
      if(pp.contains("bottom_cushion"))
        {
          pp.get("bottom_cushion", bottomCushion);
        }
      s_diffuseSolver->m_bottomSolverEpsCushion = bottomCushion;
      //      s_botSolver.m_verbosity   = 0;
      //  s_botSolver.m_numRestarts = 0;

      s_diffuseIntegratorBE  = RefCountedPtr< EBLevelBackwardEuler>( new  EBLevelBackwardEuler(grids, refRat, lev0Dom, s_diffuseOpFactory, s_diffuseSolver));
      s_diffuseIntegratorBE->setEBLG(eblgs);

      s_diffuseIntegratorTGA  = RefCountedPtr< EBLevelTGA>( new  EBLevelTGA(grids, refRat, lev0Dom, s_diffuseOpFactory, s_diffuseSolver));
      s_diffuseIntegratorTGA->setEBLG(eblgs);
    }  
}
/******************/
/******************/
void EBAMRTransport::defineFactories(bool a_atHalfTime)
{
  if (m_params.m_doDiffusion)
    {
//      MayDay::Error("diffusion set not set in EBAMRTransport");
      Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> >        >  aco(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  bco(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  bcoIrreg(nlevels);

      Vector<RefCountedPtr<EBQuadCFInterp> >                quadCFI(nlevels);

      EBAMRTransport* coarsestLevel = (EBAMRTransport*)(hierarchy[0]);
      Real           lev0Dx      = (coarsestLevel->m_dx[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblgPtr->getDomain();


      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRTransport* transportLevel = dynamic_cast<EBAMRTransport*>(hierarchy[ilev]);

          EBCellFactory fact(transportLevel->m_eblgPtr->getEBISL());
          LevelData<EBCellFAB> halfSt(transportLevel->m_eblgPtr->getDBL(),m_nComp, 4*IntVect::Unit, fact);
          if(a_atHalfTime)
            {
              transportLevel->getHalfState(halfSt);
            }
          else
            {
              Interval interv(0, m_nComp-1);
              transportLevel->m_stateOld.copyTo(interv, halfSt, interv);
            }
          transportLevel->fillCoefficients(halfSt);

          eblgs       [ilev] = *transportLevel->m_eblgPtr;
          grids       [ilev] = transportLevel->m_eblgPtr->getDBL();
          refRat      [ilev] = transportLevel->m_ref_ratio;
          aco         [ilev] = transportLevel->m_aco;
          quadCFI     [ilev] = *transportLevel->m_quadCFIPtr;
          bco         [ilev] = transportLevel->m_bco;
          bcoIrreg[ilev] = transportLevel->m_bcoIrreg;
        }

      //alpha, beta get replaced in tga solves
      Real alpha = 1;
      Real beta  = 1;
      IntVect  giv    = 4*IntVect::Unit;

      // Diffusion operator.
      int relaxType = 1;

// set the domain and the EB diffusion BCs
// for now 0 neuman
      setBCs();

      s_diffuseOpFactory =
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
         (new EBConductivityOpFactory(eblgs, quadCFI, alpha, beta, aco,
                                      bco, bcoIrreg, lev0Dx,
                                      refRat, m_domBC,
                                      m_ebBC, giv, giv, relaxType)));

    }
  (*m_bco )    .exchange(Interval(0,0));
}
/******************/
/******************/
void EBAMRTransport::
setBCs()
{
  // no diffusion of species at the domain and EB
  NeumannConductivityEBBCFactory* ebBCFact =  new NeumannConductivityEBBCFactory();
  ebBCFact->setValue(0.);
  m_ebBC = RefCountedPtr<BaseEBBCFactory>(ebBCFact);

  NeumannConductivityDomainBCFactory* domBCFact = new NeumannConductivityDomainBCFactory();
  domBCFact->setValue(0.);
  m_domBC = RefCountedPtr<BaseDomainBCFactory>(domBCFact); 
}
/******************/
/******************/
void EBAMRTransport::
fillCoefficients(const LevelData<EBCellFAB>& a_state)
{
  // for now only constant coefficients 
  CH_assert(!m_params.m_variableCoeff);
 
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      // All of these coefficients are time-independent and can be set here.
      (*m_bco)    [dit()].setVal(m_params.m_diffusionCoeff);
      (*m_bcoIrreg)[dit()].setVal(m_params.m_diffusionCoeff);
    }

  if(m_level == 0)
    {
      if((!m_params.m_variableCoeff))
        {
          pout() << "using constant diffusion coefficients"  ;
          pout() << ": diffusion coefficient = " << m_params.m_diffusionCoeff;
        }
    }
}
/******************/
/******************/
void EBAMRTransport::initialGrid(const Vector<Box>& a_new_grids)
{
  m_level_grids = a_new_grids;

  if (m_externalDriver)
    {
      CH_assert(m_isEBLGSet == true);
      CH_assert(m_eblgPtr != NULL);
      levelSetup(); 
    }
  else 
    {
      pout() << "EBAMRTransport::initialGrid not setup for self driving case" << endl;
      MayDay::Error();
    }  
}
/******************/
/******************/
void EBAMRTransport::syncWithFineLevel()
{
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if(m_hasFiner)
    {
      EBAMRTransport* finePtr = getFinerLevel();
      int nRefFine = refRatio();
      const EBLevelGrid* finer_eblgPtr = finePtr->m_eblgPtr;
      const DisjointBoxLayout& finer_dbl = finer_eblgPtr->getDBL();
      const EBISLayout& finer_ebisl      = finer_eblgPtr->getEBISL();
      // maintain flux registers
/*
      m_divFFluxRegister.define(finer_dbl,
                                m_eblgPtr->getDBL(),
                                finer_ebisl,
                                m_eblgPtr->getEBISL(),
                                m_eblgPtr->getDomain().domainBox(),
                                nRefFine,
                                m_nComp, Chombo_EBIS::instance(), s_noEBCF);

      m_scalarFluxRegister.define(finer_dbl,
                                  m_eblgPtr->getDBL(),
                                  finer_ebisl,
                                  m_eblgPtr->getEBISL(),
                                  m_eblgPtr->getDomain().domainBox(),
                                  nRefFine,
                                  m_nComp, Chombo_EBIS::instance(), s_noEBCF);
*/

// Sid added due to parallel issues: 10/12
  m_divFFluxRegister.define(*finer_eblgPtr, *m_eblgPtr,
                            nRefFine, m_nComp, s_noEBCF);

  m_scalarFluxRegister.define(*finer_eblgPtr, *m_eblgPtr,
                            nRefFine, m_nComp, s_noEBCF);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if (!s_noEBCF)
        {
          m_ebCoarToFineRedist.define(*finer_eblgPtr, *m_eblgPtr, nRefFine , m_nComp, 1);
          //define coarse to coarse redistribution object
          m_ebCoarToCoarRedist.define(*finer_eblgPtr, *m_eblgPtr, nRefFine , m_nComp, 1);
        }

      //set all the registers to zero
      if(!s_noEBCF)
        {
          m_ebCoarToFineRedist.setToZero();
          m_ebCoarToCoarRedist.setToZero();
        }
      m_divFFluxRegister.setToZero();
      m_scalarFluxRegister.setToZero();
    }
}
/******************/
/******************/
void EBAMRTransport::levelSetup()
{
  CH_assert(m_isQuadCFISet);
  CH_assert(m_isEBLGSet);

  m_nGhost = m_eblgPtr->getGhost();
  
  EBCellFactory factoryNew(m_eblgPtr->getEBISL());
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblgPtr->getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblgPtr->getDBL(),m_nComp, ivGhost, factoryNew);
  m_normalVelNew.define(m_eblgPtr->getDBL(),SpaceDim, ivGhost, factoryNew);
  m_normalVelOld.define(m_eblgPtr->getDBL(),SpaceDim, ivGhost, factoryNew);
  m_redisRHS.define(m_eblgPtr->getDBL(),m_nComp, ivGhost, factoryNew);
  EBLevelDataOps::setVal(m_redisRHS, 0.0);

  m_sets.define(m_eblgPtr->getDBL());
  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_sets[dit()] = m_eblgPtr->getEBISL()[dit()].getIrregIVS(m_eblgPtr->getDBL().get(dit()));
    }
  EBCellFactory       cellFact(m_eblgPtr->getEBISL());
  EBFluxFactory       fluxFact(m_eblgPtr->getEBISL());
  BaseIVFactory<Real> bivfFact(m_eblgPtr->getEBISL(), m_sets);

  
//  define diffusion coefficients
  // allocate A coefficients for conduction operator
  m_aco    = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, cellFact));

  // allocate species diffusion coefficient on regular and irregular cells
  m_bco     = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fluxFact));
  m_bcoIrreg= RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, bivfFact));

  EBAMRTransport* coarPtr = getCoarserLevel();
  EBAMRTransport* finePtr = getFinerLevel();       

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);

  //define redistribution object for this level
  //for now set to volume weighting
  m_ebLevelRedist.define(m_eblgPtr->getDBL(),
                         m_eblgPtr->getEBISL(),
                         m_eblgPtr->getDomain(),
                         m_nComp);

   if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      const EBLevelGrid* coEBLGPtr = coarPtr->m_eblgPtr;

      m_ebCoarseAverage.define(m_eblgPtr->getDBL(),
                               coEBLGPtr->getDBL(),
                               m_eblgPtr->getEBISL(),
                               coEBLGPtr->getEBISL(),
                               coEBLGPtr->getDomain().domainBox(),
                               nRefCrse,
                               m_nComp, Chombo_EBIS::instance());
      m_ebFineInterp.define(m_eblgPtr->getDBL(),
                            coEBLGPtr->getDBL(),
                            m_eblgPtr->getEBISL(),
                            coEBLGPtr->getEBISL(),
                            coEBLGPtr->getDomain().domainBox(),
                            nRefCrse,
                            m_nComp);

      // maintain levelgodunov
      m_ebLevelGodunov.define(m_eblgPtr->getDBL(),
                              coEBLGPtr->getDBL(),
                              m_eblgPtr->getEBISL(),
                              coEBLGPtr->getEBISL(),
                              m_eblgPtr->getDomain(),
                              nRefCrse,
                              m_dx,
                              m_hasCoarser,
                              m_hasFiner,
                              m_ebPatchGodunov);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if(!s_noEBCF)
        {
          CH_TIME("fineToCoar_defs");
/*
          m_ebFineToCoarRedist.define(*m_eblgPtr, *coEBLGPtr,
                                      nRefCrse, m_nComp,
                                      m_params.m_redistRad);
*/
// new define Sid: 10/12
          m_ebFineToCoarRedist.define(m_eblgPtr->getDBL(), 
                                      coEBLGPtr->getDBL(),
                                      m_eblgPtr->getEBISL(),
                                      coEBLGPtr->getEBISL(),
                                      coEBLGPtr->getDomain().domainBox(),
                                      nRefCrse,
                                      m_nComp, m_params.m_redistRad);          
          m_ebFineToCoarRedist.setToZero();
        }

      int nvarQuad = 1; 

/*
      m_quadCFI = RefCountedPtr<EBQuadCFInterp>
        (new EBQuadCFInterp(m_eblgPtr->getDBL(),
                            coEBLGPtr->getDBL(),
                            m_eblgPtr->getEBISL(),
                            coEBLGPtr->getEBISL(),
                            coEBLGPtr->getDomain(),
                            nRefCrse, nvarQuad,
                            (*m_eblgPtr->getCFIVS()),
                            Chombo_EBIS::instance()));
*/
      coarPtr->syncWithFineLevel();
    }
  else
    {
//      m_quadCFI = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp());
      m_ebLevelGodunov.define(m_eblgPtr->getDBL(),
                              DisjointBoxLayout(),
                              m_eblgPtr->getEBISL(),
                              EBISLayout(),
                              m_eblgPtr->getDomain(),
                              m_ref_ratio,
                              m_dx,
                              m_hasCoarser,
                              m_hasFiner,
                              m_ebPatchGodunov);
    }

  m_sets.define(m_eblgPtr->getDBL());
  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_eblgPtr->getDBL().get(dit());
      m_sets[dit()] = m_eblgPtr->getEBISL()[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_eblgPtr->getEBISL(), m_sets);
  //the ghost cells on the mass redistribution array
  //are tied to the redistribution radius
  m_massDiff.define(m_eblgPtr->getDBL(), m_nComp, ivGhost, factory);
  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
    }
}
/******************/
/******************/
EBAMRTransport*
EBAMRTransport::
getCoarserLevel() const
{
  EBAMRTransport* retval = NULL;
  if(m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRTransport*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/******************/
/******************/
EBAMRTransport*
EBAMRTransport::
getFinerLevel() const
{
  EBAMRTransport* retval = NULL;
  if(m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast<EBAMRTransport*>(m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/******************/ 
/******************/
void EBAMRTransport::initialData()
{

  const EBPhysIBC* const ebphysIBCPtr = m_ebPatchGodunov->getEBPhysIBC();

  ebphysIBCPtr->initialize(m_stateNew, m_eblgPtr->getEBISL());
  ebphysIBCPtr->initialize(m_stateOld, m_eblgPtr->getEBISL());

}
/******************/
/******************/
void EBAMRTransport::postInitialize()
{
    // Average the finer grid levels to be consistent with this one.
  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRTransport* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }

/*
  // only if the equation being solved conserves something
  // Here we record the totals for the conserved quantities on grid level 0.
  if(!m_hasCoarser)
    {
      // Sum the quantities.
      int densityIndex = m_ebPatchGodunov->densityIndex();
      sumConserved(m_originalMass, densityIndex);
    }
//    if(m_level==0) defineSolvers(); 
*/
}
/******************/
/******************/
void EBAMRTransport::
sumConserved(Real& a_sumcons,
             const int& a_ivar) const
{
  Real sumgrid = 0;
  for(DataIterator dit= m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& state = m_stateNew[dit()];
      Box thisBox = m_eblgPtr->getDBL().get(dit());
      IntVectSet uberIVS(thisBox);
      const EBISBox& ebisBox = m_eblgPtr->getEBISL()[dit()];
      for(VoFIterator vofit(uberIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real consVal = state(vof, a_ivar);
          Real volFrac = ebisBox.volFrac(vof);
          Real volume = volFrac;
          sumgrid += consVal*volume;
        }
    }

  Vector<Real> all_sum;
  gather(all_sum,sumgrid,uniqueProc(SerialTask::compute));
  Real sumallgrid = 0.;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < all_sum.size(); ++i)
        {
          sumallgrid += all_sum[i];
        }
    }
  broadcast(sumallgrid,uniqueProc(SerialTask::compute));
  a_sumcons = sumallgrid;
}
/******************/
/******************/
Real EBAMRTransport::computeDt()
{
  // dt is provided by external driver
  m_dt = *m_dtPtr;
  return m_dt;
}
/******************/
/******************/
Real EBAMRTransport::computeInitialDt()
{
  // dt is provided by external driver
  m_dt = *m_dtPtr;
  return m_dt;
}
/******************/
/******************/
void EBAMRTransport::assignAdvVelPtr(const LevelData<EBFluxFAB>* a_advVelPtr,
                                     const LayoutData< Vector< BaseIVFAB<Real> * > >* a_coveredAdvVelLoPtr,
                                     const LayoutData< Vector< BaseIVFAB<Real> * > >* a_coveredAdvVelHiPtr)
{
  m_advVelPtr          = a_advVelPtr;
  m_coveredAdvVelLoPtr = a_coveredAdvVelLoPtr;
  m_coveredAdvVelHiPtr = a_coveredAdvVelHiPtr;
}
/******************/
/******************/
void EBAMRTransport::
setNormalVelOld(const LevelData<EBCellFAB>& a_normalVel)
{
  a_normalVel.copyTo(m_normalVelOld);
}
/******************/
/******************/
void EBAMRTransport::
setNormalVelOld(const LevelData<EBFluxFAB>& a_advVel,
                const LayoutData< Vector< BaseIVFAB<Real> * > >& a_coveredAdvVelLo,
                const LayoutData< Vector< BaseIVFAB<Real> * > >& a_coveredAdvVelHi)
{
  m_oldNormalVelSet = true;
  m_ebLevelGodunov.computeNormalVel(m_normalVelOld, a_advVel, a_coveredAdvVelLo, a_coveredAdvVelHi);
}
/******************/
/******************/
#ifdef CH_USE_HDF5
/******************/
/******************/
void EBAMRTransport::writePlotHeaderOld    (HDF5Handle& a_handle) const
{
}
/******************/
/******************/
void EBAMRTransport::writePlotLevelOld     (HDF5Handle& a_handle) const
{
}
/******************/
/******************/
void EBAMRTransport::writePlotHeader       (HDF5Handle& a_handle) const
{
}
/******************/
/******************/
void EBAMRTransport::writePlotLevel        (HDF5Handle& a_handle) const
{
}
/******************/
/******************/
void EBAMRTransport::writeCheckpointHeader (HDF5Handle& a_handle) const
{
  //stuff in checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman 
}
/******************/
/******************/
void EBAMRTransport::writeCheckpointLevel  (HDF5Handle& a_handle) const
{
 
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);
 
  // write data for this level
  write(a_handle, m_stateOld, "dataOld");
  write(a_handle, m_stateNew, "dataNew");
}
/******************/
/******************/
void EBAMRTransport::readCheckpointHeader  (HDF5Handle& a_handle)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRTransport::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman. 
}
/******************/
/******************/
void EBAMRTransport::readCheckpointLevel   (HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  // only the case driven by external driver for now
  CH_assert(m_externalDriver);
  CH_assert(m_isEBLGSet);
  CH_assert(m_eblgPtr != NULL);

  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRTransport::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

// refRat is alread set in define
/*
  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];
*/

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

/*
  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];
*/

// already checked in the driver class
/*
  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids);
  if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }
*/

  //this keeps the order of the AMRLevel m_level_grids
  //consistent with m_eblg.getDBL()
  LayoutIterator lit = m_eblgPtr->getDBL().layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = m_eblgPtr->getDBL().get(lit());
      m_level_grids.push_back(b);
    }

// this is in levelSetup
/*
  EBCellFactory factoryNew(m_eblgPtr->getEBISL());
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblgPtr->getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblgPtr->getDBL(),m_nComp, ivGhost, factoryNew);
  m_normalVelNew.define(m_eblgPtr->getDBL(),SpaceDim, ivGhost, factoryNew);
  m_normalVelOld.define(m_eblgPtr->getDBL(),SpaceDim, ivGhost, factoryNew);
  m_redisRHS.define(m_eblgPtr->getDBL(),m_nComp, ivGhost, factoryNew);
  EBLevelDataOps::setVal(m_redisRHS, 0.0);
*/
  // Set up data structures
  levelSetup();
//  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNew",
                                      m_eblgPtr->getDBL(),
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOld",
                                      m_eblgPtr->getDBL(),
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }

}
/******************/
/******************/
#endif

#include "NamespaceFooter.H"

