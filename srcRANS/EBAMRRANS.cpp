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

#include "EBAMRRANS.H"
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
#include "EBPatchGodunov.H"

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

Vector<RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > > EBAMRRANS::s_diffuseOpFactory;

Vector<RefCountedPtr< EBLevelBackwardEuler> >                     EBAMRRANS::s_diffuseIntegratorBE;

Vector<RefCountedPtr< EBLevelTGA> >                               EBAMRRANS::s_diffuseIntegratorTGA;

Vector<RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > > > EBAMRRANS::s_diffuseSolver;

//BiCGStabSolver<LevelData<EBCellFAB> >                  EBAMRRANS::s_botSolver;

EBSimpleSolver                 EBAMRRANS::s_botSolver;
bool EBAMRRANS::s_solversDefined = false;
/******************/
EBAMRRANS::
~EBAMRRANS()
{
  m_ebPatchModel = RefCountedPtr<EBPatchRANSModel>();

  if ((m_level == 0) && m_params.m_doDiffusion)
   {
     clearSolvers();
   }
}
/******************/
EBAMRRANS::clearSolvers()
{
  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
  {
    s_diffuseOpFactory[iEqn] = RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
    s_diffuseIntegratorBE[iEqn] = RefCountedPtr< EBLevelBackwardEuler>();
    s_diffuseIntegratorTGA[iEqn] = RefCountedPtr< EBLevelTGA>();
    s_diffuseSolver[iEqn] = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >();
  }  
}
/******************/
EBAMRRANS::
EBAMRRANS(const EBAMRRANSParams& a_params,
          const RefCountedPtr<EBPatchRANSModelFactory>& a_modelFactory,
          bool  a_externalDriver):
    m_params(a_params),
    m_ebPatchRANSModelFactory(a_modelFactory),
    m_externalDriver(a_externalDriver)
{
  m_isDxSet = false;
  m_isDBLSet = false;
  m_isEBLGSet = false;
  m_isQuadCFISet = false;
  m_isEBISLSet = false;
  m_eblgPtr = NULL;
  m_oldNormalVelSet = false; 
}
/******************/
EBAMRRANS::
define(EBAMRLevel*          a_coarser_level_ptr,
       const ProblemDomain& a_problem_domain,
       int                  a_level,
       int                  a_ref_ratio)
{
  if (m_params.m_verbosity >=3)
    {
      pout() << "EBAMRRANS::define, level=" << a_level << endl;
    }

  EBPatchGodunov::useConservativeSource(true);
  m_isDefined = true;
  EBAMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);

  if (a_coarser_level_ptr != NULL)
    {
      EBAMRRANS* amrg_ptr =
        dynamic_cast<EBAMRRANS*>(a_coarser_level_ptr);
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
      pout() << "EBAMRRANS class is not setup to drive itself" << endl;
      MayDay::Error();
    }

    m_nGhost = 4; // hmm...

  m_ebPatchModel = RefCountedPtr<EBPatchRANSModel>();
  m_ebPatchModel = RefCountedPtr<EBPatchRANSModel>(m_ebPatchModelFactory->create());
  m_ebPatchModel->define(m_problem_domain, m_dx);

  m_nComp      = m_ebPatchModel->numConserved();
  m_nPrim      = m_ebPatchModel->numPrimitives();
  m_stateNames = m_ebPatchModel->stateNames();
  m_primNames  = m_ebPatchModel->primNames();
  m_ref_ratio  = a_ref_ratio;
}
/**********/
void EBAMRRANS::assignDx(RealVect a_dx)
{
  m_isDxSet = true;
  m_dx = a_dx;
}
/**********/
Real EBAMRRANS::advance()
{
  pout() << "advancing EBAMRRANS with dt = " << m_dt << endl;

  EBPatchGodunov::s_whichLev = m_level;

  if(m_params.m_variableCoeff && (m_level==0) && m_params.m_doDiffusion) defineSolvers();

  m_stateNew.copyTo(m_stateNew.interval(), m_stateOld, m_stateOld.interval());

  EBCellFactory fact(m_eblgPtr->getEBISL());
  LevelData<EBCellFAB>      divergeF(m_eblgPtr->getDBL(),m_nComp          , 4*IntVect::Unit, fact);

  pout() << "EBAMRRANS::getting diverence of flux" << endl;

  fluxDivergence( divergeF);

  if(!m_params.m_doDiffusion)
    {
      explicitAdvance(divergeF);
    }

  else
    {
//      pout() << "EBAMRRANS::diffusion is not set" << endl;
      // diffusion dance
      pout() << "EBAMRRANS::step1: getting Ustar (U^n + explicit contribution)" << endl;

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

      // add production-dissipation
      pout() << "step 2(i): adding producion - dissipation to transport equation"
      LevelData<EBCellFAB> netSource(m_eblg->getDBL(), m_nEqn, 4*IntVect::Unit, fact);
      addNetSource(netSource, UStar); 
    
      // now add the div(kappaGradS)
      pout() << "step2: solve for diffusion term and add to transport equation"
      LevelData<EBCellFAB> dtDivDGradS(m_eblgPtr->getDBL(),       1, 4*IntVect::Unit, fact);
      getDivDGradS(dtDivDGradS, UStar);

      pout() << "putting state into m_statenew and flooring" << endl;
      finalAdvance(UStar);
    }

  if(m_params.m_checkMaxMin)
    {
      CH_TIME("max_min_check");
      pout() << "EBAMRRANS::after advance max mins for data for level " << m_level << endl;
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
/**********/
void EBAMRRANS::fluxDivergence(LevelData<EBCellFAB>& a_divergeF)
{
  CH_assert(a_divergeF.nComp() == m_nComp);
  CH_assert(m_externalDriver);
  CH_assert(m_oldNormalVelSet);
  CH_assert(m_oldTimeSet);
  CH_assert(m_advVelPtr != NULL);
  CH_assert(m_coveredAdvVelLoPtr != NULL);
  CH_assert(m_coveredAdvVelHiPtr != NULL);

  pout() << "EBAMRRANS::taking flux divergence on level " << m_level << " with dt = " << m_dt << endl;
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

  m_ebLevelRANS.computeNormalVel(m_normalVelNew, *m_advVelPtr, *m_coveredAdvVelLoPtr, *m_coveredAdvVelHiPtr);

  if(m_hasCoarser)
    {
      EBAMRRANS* coarPtr = getCoarserLevel();
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

  m_ebLevelRANS.divergeF(a_divergeF,
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
/**********/
void EBAMRRANS::
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
          EBCellFactory fact(m_eblgPtr->getEBISL());
          LevelData<EBCellFAB>  hypSrc(m_eblgPtr->getDBL(), m_nComp,4*IntVect::Unit, fact);
          explicitHyperbolicSource(hypSrc, m_stateOld, true);

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
/**********/
void EBAMRRANS::
explicitHyperbolicSource(LevelData<EBCellFAB>&       a_hypSource,
                         const LevelData<EBCellFAB>& a_state,
                         bool a_doNormalization)
{
  pout() << "EBAMRRANS::computing explicit hyperbolic source" << endl;
  EBCellFactory fact(m_eblgPtr->getEBISL());

  LevelData<EBCellFAB>* stateCoar = NULL;
  getCoarserState(stateCoar);

  LevelData<EBCellFAB> netSource(m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);
  LevelData<EBCellFAB> diffSource(m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);

  kappaDiffusiveSource(diffSource, a_state, stateCoar);
  kappaNetSource(netSource, a_state, stateCoar); // production - dissipation
  
  EBLevelDataOps::sum(a_hypSource,diffSource,netSource);

  if (m_hasCoarser)
    {
      delete stateCoar;
    }

  if (a_doNormalization)
    {
      //finally all these things are multiplied by kappa so we need to make 
      //them normalized 
      KappaSquareNormal normalizinator(*m_eblgPtr);
      normalizinator(a_hypSource);
    }
}
/**********/
void EBAMRRANS::
getCoarserState(LevelData<EBCellFAB>* & a_stateCoar)
{
  LevelData<EBCellFAB>* stateCoar = NULL;
  if(m_hasCoarser)
    {
      EBAMRRANS* coarPtr = getCoarserLevel();
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
/**********/
EBAMRRANS* EBAMRRANS::getCoarserLevel() const
{
  EBAMRRANS* retval = NULL;
  if(m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRRANS*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/**********/
EBAMRRANS* EBAMRRANS::getFinerLevel() const
{
  EBAMRRANS* retval = NULL;
  if(m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast<EBAMRRANS*>(m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/*********/
/// set  output to (del dot (kappa grad S))  
void
EBAMRRANS::
kappaDiffusiveSource(LevelData<EBCellFAB>& a_kappaDiffSource,
                     const LevelData<EBCellFAB>& a_state,
                     const LevelData<EBCellFAB>* a_stateCoar)
{
  //set the a coefficients to time n

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblgPtr->getDBL().get(dit());
      for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
        {
          (*m_aco[iEqn])[dit()].setVal(m_params.m_alphaCoeff); // caution!!! hardwiring to be same constant a coeff for all eqns
        }
    }

  // Compute the diffusion term = volfrac(div kappa grad S)
  Real alpha = 0;   Real beta = 1; //want just the div(flux) part of the operator
  int nghost = 4;
  EBCellFactory fact(m_eblgPtr->getEBISL());

  // eqation after eqation
  for (int iEqn = 0; iEqn < m_nEqn; iEqn++) 
   {
     LevelData<EBCellFAB> kappaSrc(m_eblgPtr->getDBL(),  1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB> eqn(m_eblgPtr->getDBL(), 1, nghost*IntVect::Unit, fact);
     LevelData<EBCellFAB>* eqnCoar = NULL;
     Interval eqnInt(iEqn, iEqn);
     Interval thisInt(0, 0);
     a_state.copyTo(eqnInt, eqn, thisInt);

     if(m_hasCoarser)
      {
        eqnCoar = new LevelData<EBCellFAB>(m_eblgPtr->getDBL(),  1, nghost*IntVect::Unit, fact);
        a_stateCoar->copyTo(eqnInt, *eqnCoar, thisInt);
      }

     s_diffuseIntegratorBE[iEqn]->applyOperator(kappaSrc, eqn, eqnCoar, m_level, alpha, beta, true);

     if (m_hasCoarser)
      {
        delete eqnCoar
      }

     kappaSrc.copyTo(thisInt, a_kappaDiffSource, eqnInt);
   } // end for iEqn
}
/*********/
void EBAMRRANS::
kappaNetSource(LevelData<EBCellFAB>& a_kappaNetSource,
               const LevelData<EBCellFAB>& a_state,
               const LevelData<EBCellFAB>* a_stateCoar)
{
  m_ebLevelRANS.computeNetSource(a_kappaNetSource, a_state, m_timeOld, m_dtOld);
  EBLevelDataOps::kappaWeight(a_kappaNetSource);
}
/*********/
void EBAMRRANS::
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
/*********/
void EBAMRRANS::
explicitAdance(const LevelData<EBCellFAB>& a_divergeF)
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

  m_ebLevelRANS.floorPrimitives(m_stateNew, m_timeOld, m_dtOld);
}
/*********/
void EBAMRRANS::
hyperbolicRedistribution(LevelData<EBCellFAB>& a_state)
{
  //explicit redistribution here because these are the hyperbolic terms

  m_ebLevelRedist.setToZero();
/*
  if(m_params.m_useMassRedist)
    {
      //if use mass weighting, need to
      //fix weights of redistribution object
      int densityIndex = m_ebPatchModel->densityIndex();
      a_state.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(a_state, densityIndex);
    }
*/
  Interval consInterv(0, m_nComp-1);
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), consInterv);
    }
  //

  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
    {
      if(m_params.m_useMassRedist)
        {
          //if use mass weighting, need to
          //fix weights of redistribution object
//          int densityIndex = m_ebPatchModel->densityIndex();
          a_state.exchange(Interval(iEqn,iEqn));
          m_ebLevelRedist.resetWeights(a_state, iEqn);
        }

      m_ebLevelRedist.redistribute(a_state, Interval(iEqn,iEqn));
    }
  m_ebLevelRedist.setToZero();
}
/*********/
//U* = Un + dt*L^c(U^n) + (dt/2)*L^v(U^n)
void EBAMRRANS::
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

  m_ebLevelRANS.floorConserved(a_UStar, m_time, m_dt);
}
/*********/
void EBAMRRANS::
addNetSource(LevelData<EBCellFAB>& a_netSource,
             LevelData<EBCellFAB>& a_UStar)
{
  LevelData<EBCellFAB>* stateCoar = NULL;
  getCoarserState(stateCoar);
  EBCellFactory fact(m_eblgPtr->getEBISL());
  int nghost = 4;
  LevelData<EBCellFAB> kappaConsNetSource(m_eblgPtr->getDBL(), m_nComp, nghost*IntVect::Unit, fact);
  LevelData<EBCellFAB> nonConsNetSource(m_eblgPtr->getDBL(), m_nComp, nghost*IntVect::Unit, fact);

  kappaNetSource(kappaConsNetSource, a_UStar, stateCoar);

  // get the non conservative version
  EBLevelDataOps::setToZero(nonConsNetSource);
  EBLevelDataOps::incr(nonConsNetSource, kappaConsNetSource, 1.0);
  KappaSquareNormal normalizinator(*m_eblgPtr);
  normalizinator(nonConsNetSource);

  updateStateByNetSourceAndRedistribute(a_netSource, kappaConsNetSource, nonConsNetSource, a_UStar);

  if(m_hasCoarser)
    {
      delete stateCoar;
    }
}
/*********/
void EBAMRRANS::
updateStateByNetSourceAndRedistribute(LevelData<EBCellFAB>& a_netSource,
                                      LevelData<EBCellFAB>& a_kappaConsNetSource,
                                      LevelData<EBCellFAB>& a_nonConsNetSource,
                                      LevelData<EBCellFAB>& a_UStar)
{
  // make net source = kappa*cons + (1-kappa)*noncons
  EBLevelDataOps::setToZero(a_netSource);
  EBLevelDataOps::incr(a_netSource, a_kappaConsNetSource, 1.0);
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblgPtr->getDBL().get(dit());
      IntVectSet ivs = m_eblgPtr->getEBISL()[dit()].getIrregIVS(region);
      for(VoFIterator vofit(ivs, m_eblgPtr->getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Real kappa   =  m_eblg.getEBISL()[dit()].volFrac(vofit());
          for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
            {
              Real kapCons =  a_kappaConsNetSource[dit()](vofit(), iEqn);
              Real nonCons =  a_nonConsNetSource[dit()](vofit(), iEqn);
              a_netSource[dit()](vofit(), iEqn) = kapCons + (1.-kappa)*nonCons;
              m_massDiff[dit()](vofit(), iEqn) = m_dt*((1.-kappa)*kapCons - kappa*(1.-kappa)*nonCons);
            }
        }
    }
  Interval interv(0, m_nEqn-1);
  for (DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB incr(m_eblgPtr->getEBISL()[dit()], m_eblg.getDBL().get(dit()), m_nEqn);
      incr.setVal(0.);
      //this makes incr = dissipation function = divsigma u
      incr += a_netSource[dit()];
      //this makes incr = dt*(netSource)
      incr *= m_dt;
      //this adds effects of net source into the state
      a_UStar[dit()] += incr;

      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), interv);
      m_massDiff[dit()].setVal(0.);
    }

  //this smushes the source difference into the state

  a_UStar.exchange(interv);

  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
    {
      m_ebLevelRedist.resetWeights(a_Ustar, iEqn);
      m_ebLevelRedist.redistribute(m_redisRHS, Interval(iEqn,iEqn));
    }

  m_ebLevelRedist.setToZero();
}
/*********/
void EBAMRRANS::
getDivDGradS(LevelData<EBCellFAB>& a_divDGradS,
                    LevelData<EBCellFAB>& a_Ustar)
{
  EBCellFactory fact(m_eblgPtr->getEBISL());

  EBFluxRegister*       coarFRPtr = NULL;
  EBFluxRegister*       fineFRPtr = NULL;
  LevelData<EBCellFAB>* UCoarOldPtr = NULL;
  LevelData<EBCellFAB>* UCoarNewPtr = NULL;
  Real tCoarOld= 0;
  Real tCoarNew= 0;

  if(m_hasCoarser)
    {
      // Get coarser data if we have it.
      EBAMRRANS* coarRANS = getCoarserLevel();
      coarFRPtr = &coarRANS->m_scalarFluxRegister;
      tCoarNew = coarRANS->m_time;
      tCoarOld = tCoarNew - coarRANS->m_dt;

      const EBLevelGrid& ceblg = *coarRANS->m_eblgPtr;
      EBCellFactory cfact(ceblg.getEBISL());
      UCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), m_nEqn, 4*IntVect::Unit, cfact);
      UCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), m_nEqn, 4*IntVect::Unit, cfact);

      coarRANS->m_stateOld.copyTo(*UCoarOldPtr);
      coarRANS->m_stateNew.copyTo(*UCoarNewPtr);
    }

  if (m_hasFiner) 
    {

      // Get finer flux registers if we have them.

      fineFRPtr = &m_scalarFluxRegister;
    }

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblgPtr->getDBL().get(dit());
      for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
        {
          (*m_aco[iEqn])[dit()].setVal(m_params.m_alphaCoeff); // caution!!! hardwiring to be same constant a coeff for all eqns
        }
    }

  // eqation after eqation
  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
   {
     LevelData<EBCellFAB> Sold(m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fact);
     LevelData<EBCellFAB> Snew(m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fact);
     LevelData<EBCellFAB> rhs(m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fact);
     EBLevelDataOps::setToZero(rhs);
     
     LevelData<EBCellFAB>* SCoarOld = NULL;
     LevelData<EBCellFAB>* SCoarNew = NULL;
     Interval eqnInt(iEqn, iEqn);
     Interval thisInt(0, 0);
     a_UStar.copyTo(eqnInt, Sold, thisInt);
     
     if (m_hasCoarser)
       {
         UCoarOld->copyTo(eqnInt, *SCoarOld, thisInt);
         UCoarNew->copyTo(eqInt, *SCoarNew, thisInt);
       }
    
     m_redisRHS.copyTo(eqnInt, rhs, dstInt);
     EBLevelDataOps::setToZero(m_redisRHS);

     if (m_params.m_backwardEuler)
      {
        s_diffuseIntegratorBE[iEqn]->updateSoln(Snew, Sold, rhs, fineFRPtr, coarFRPtr, SCoarOldPtr, SCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt/2., m_level, true, iEqn);
      }
     else
      {
        s_diffuseIntegratorTGA[iEqn]->updateSoln(Snew, Sold, rhs, fineFRPtr, coarFRPtr, SCoarOldPtr, SCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt/2., m_level, true, iEqn);
      }

     EBLevelDataOps::scale(Sold, -1.0);
     EBLevelDataOps::scale(Snew,  1.0);

     for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
      {
        a_dtDivDGradS[dit()].setVal(iEqn, 0.);
        //sets divDgradS = ((Snew-Sold)/dt)*dt
        int isrc = 0; int idst = iEqn; int inco = 1;
        a_dtDivDGradS[dit()].plus(Sold[dit()], isrc, idst, inco); // negative sign is in the call above
        a_dtDivDGradS[dit()].plus(Snew[dit()], isrc, idst, inco); 
        
        //this makes divDgradS = (alpha (Snew-Sold)/dt)*dt
        a_dtDivDGradS[dit()].mult((*m_aco[iEqn])[dit()], isrc, idst, inco);

        a_UStar[dit()].plus(a_dtDivDGradS[dit()], idst, idst, inco);
      }

     if (m_hasCoarser)
      {
        delete SCoarOldPtr;
        delete SCoarNewPtr;
      }
   } // end for iEqn  
  
  if (m_hasCoarser)
    {
      delete UCoarOldPtr;
      delete UCoarNewPtr,
      delete coarFRPtr;
    }

  if (m_hasFiner)
    {
      delete fineFRPtr;
    }
}
/**********/
void EBAMRRANS::
finalAdvance(LevelData<EBCellFAB>& a_Ustar)
{
  //set stateNew to udbst
  Interval comps(0, m_nComp-1);
  a_Ustar.copyTo(comps, m_stateNew, comps);

  m_ebLevelRANS.floorConserved(m_stateNew, m_time, m_dt);
}
/**********/
void EBAMRRANS::
postTimeStep()
{
//  pout() << " in EBAMRRANS postTimeStep for level " << m_level << endl;

 // Average the finer grid levels to be consistent with this one.

  m_oldNormalVelSet = false;

  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRRANS* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }

 // this does the refluxing and redistribution evil dance
  postTimeStepRefluxRedistDance();
  m_ebLevelRANS.floorPrimitives(m_stateNew, m_time, m_dt);

  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
      m_redisRHS[dit()].setVal(0.0);
    }
}
/*********/
void EBAMRRANS::postTimeStepRefluxRedistDance()
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
/*********/
void EBAMRRANS::
resetWeights()
{
  // currently resetting weight based on only one equation 
  // I think its reasonable
  if (m_params.m_useMassRedist)
    {
      int weighingvar = m_ebPatchModel->stencilWeightEqnIndex();
      m_stateNew.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(m_stateNew, weighingvar);

      if(m_hasCoarser && (!s_noEBCF))
        {
          EBAMRRANS* coarPtr = getCoarserLevel();
          coarPtr->m_stateNew.exchange(Interval(0, m_nComp-1));
          m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, weighingvar);
        }
      if(m_hasFiner && (!s_noEBCF))
        {
          m_ebCoarToFineRedist.resetWeights(m_stateNew, weighingvar);
          m_ebCoarToCoarRedist.resetWeights(m_stateNew, weighingvar);
        }
    }
}
/*********/
void EBAMRRANS::refluxRedistInteraction()
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
/**********/
void EBAMRRANS::coarseFineRedistribution(const Interval& a_interv)
{
  if(m_hasCoarser)
    {
      if(m_params.m_doSmushing && (!s_noEBCF))
        {
          //redistibute to coarser level
          EBAMRRANS* coarPtr = getCoarserLevel();
          //put mass directly into state
          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, a_interv);
          m_ebFineToCoarRedist.setToZero();
        }
    }
  if (m_hasFiner)
    {
      EBAMRRANS* finePtr = getFinerLevel();
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
/*********/
void EBAMRRANS::explicitReflux(const Interval& a_interv)
{
  if (m_hasFiner)
    {
      Real scale = -1.0/m_dx[0];
      m_divFFluxRegister.reflux(m_stateNew, a_interv, scale);
      m_divFFluxRegister.setToZero();
    }
}
/*********/
void EBAMRRANS::refluxRHS()
{
  //this does the refluxing and redistribution evil dance
  Interval interv(0, m_nComp-1);
  if (m_hasFiner)
    {
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];
      Interval interflux(0, m_nComp-1);
      m_scalarFluxRegister.reflux(m_redisRHS, interv, interflux, scale);

      m_scalarFluxRegister.setToZero();
    }
}
/*********/
void EBAMRRANS::implicitReflux()
{
  pout() << "EBAMRRANS::implicit redist/reflux for level" << m_level << endl;
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
          EBAMRRANS* RANSLevel = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);

          EBCellFactory fact(RANSLevel->m_eblgPtr->getEBISL());
          deltaScalar[ilev]    = new LevelData<EBCellFAB>(RANSLevel->m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);
          dtRefluxDivergeS[ilev] = new LevelData<EBCellFAB>(RANSLevel->m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, fact);

          //copy redistRHS to reflux divergence holders
          Interval srcComp, dstComp;
          srcComp = Interval(0, m_nComp-1);
          dstComp = Interval(0, m_ncomp-1);
          RANSLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeS[ilev], dstComp);

          EBLevelDataOps::scale(*dtRefluxDivergeS[ilev], RANSLevel->m_dt);

          //these MATTER especially with subcycling
          EBLevelDataOps::setVal(*deltaScalar[ilev], 0.0);
        }

      pout() << "getting scalar reflux/redist increment"  << endl;
      int baseLev = Max(m_level-1, 0);
      Real baseDt = ((EBAMRRANS*)hierarchy[baseLev])->m_dt;

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
          EBAMRRANS* RANSLevel = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);
          EBLevelDataOps::setToZero(RANSLevel->m_redisRHS);
        }
    }
}
/*********/
int EBAMRRANS::
getFinestLevel()
{
  int imax = 0;
  Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
  for(int ilev = 0; ilev < hierarchy.size(); ilev++)
    {
      EBAMRRANS* RANSLevel = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);
      int numBoxes = 0;
      if(RANSLevel->m_eblgPtr->getDBL().isClosed()) numBoxes = RANSLevel->m_eblgPtr->getDBL().size();
      if(numBoxes > 0)
        imax++;
      else
        break;
    }
  return imax-1;
}
/*********/
void EBAMRRANS::
getRefluxDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeS,
                int a_baseLev, int a_finestLev, Real a_dtBase)
{
  Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for(int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRRANS* RANSLevel = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);
      Vector<RefCountedPtr<LevelData<EBCellFAB> > >& aco = RANSLevel->m_aco;
       
      for(DataIterator dit = RANSLevel->m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = RANSLevel->m_eblgPtr->getDBL().get(dit());
          for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
            {
              (*aco[iEqn])[dit()].setVal(m_params.m_alphaCoeff); // caution!!! hardwiring to be same constant a coeff for all eqns
            }
        }
    }

  // eqation after eqation
  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
   {
     Interval eqnInt(iEqn, iEqn);
     Interval thisInt(0, 0);

     //set alpha to 1, beta = -dtbase
     s_diffuseIntegratorBE[iEqn]->resetSolverAlphaAndBeta(1.0, -a_dtBase);
 
     //the sad way of doing it :-(
     Vector<LevelData<EBCellFAB>* > thisDeltaS(a_finestLev+1, NULL);
     Vector<LevelData<EBCellFAB>* > thisDtRefluxDivergeS(a_finestLev+1, NULL);

     for (int ilev = 0; ilev <= a_finestLev; ilev++)
      {
        EBAMRRANS* RANSLevel = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);
        EBCellFactory fact(RANSLevel->m_eblgPtr->getEBISL());
        thisDeltaS[ilev]        = new LevelData<EBCellFAB>(RANSLevel->m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fact);  
        thisDtRefluxDivergeS[ilev] = new LevelData<EBCellFAB>(RANSLevel->m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fact);
    
        a_dtRefluxDivergeS[ilev]->copyTo(eqnInt, *thisDtRefluxDivergeS, thisInt);
        a_deltaScalar[ilev]->copyTo(eqnInt, *thisDeltaS, thisInt); 
      }

     //solve equation (rho C_v I - dt Lv) delta = dt*Dr(Fm)
     //rhs already multiplied by dt
     //first true = zero phi
     //second true = force homogeneous bcs (solving for a delta t)
     s_diffuseSolver[iEqn]->solve(thisDeltaS, thisDtRefluxDivergeS, a_finestLev, a_baseLev, true, true);

     for (int ilev = 0; ilev <= a_finestLev; ilev++)
      {
        thisDeltaS[ilev]->copyTo(thisInt, *a_deltaScalar[ilev], eqnInt);

        delete thisDeltaS[ilev];
        delete thisDtRefluxDivergeS[ilev];
      }

   } // end for iEqn
}
/*********/
// scalar += dt*Ls(deltaS)
void EBAMRRANS::
incrScalarByDeltaS(Vector<LevelData<EBCellFAB>* >& a_deltaScalar,
                    int a_baseLev, int a_finestLev)
{
  Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRRANS* RANSLevel   = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);
      RANSLevel->incrScalarByDeltaS(RANSLevel->m_stateNew, *a_deltaScalar[ilev]);
    }
}
/*********/ 
// scalar += dt*Ls(deltaS)
void EBAMRRANS::
incrScalarByDeltaS(LevelData<EBCellFAB>& a_state,
                   LevelData<EBCellFAB>& a_deltaScalar)
{
  for(DataIterator dit = m_eblgPtr->getDBL().dataIterator(); dit.ok(); ++dit)
    {
      //add change in scalar to the state
      {
        int isrc = 0; int idst = 0; int inco = m_nEqn;
        a_state[dit()].plus(a_deltaScalar[dit()], isrc, idst, inco);
      }
    }
}
/*********/
void EBAMRRANS:: 
tagCells(IntVectSet& a_tags)
{
}
/********/
void EBAMRRANS::
tagCellsInit(IntVectSet& a_tags)
{
}
/********/
void EBAMRRANS::
regrid(const Vector<Box>& a_new_grids)
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
  //      eblgNew is set through the pointer: not safe, use caution

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
      EBAMRRANS* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew,
                                 coarPtr->m_stateNew,
                                 interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
  m_stateNew.copyTo(interv,m_stateOld, interv);
}
/*********/
void EBAMRRANS::levelSetup()
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

  // Allocations
  m_aco.resize(m_nEqn);
  m_bco.resize(m_nEqn);
  m_bcoIrreg.resize(m_nEqn);

//  define diffusion coefficients 
  for (int iEqn = 0; iEqn <  m_nEqn; iEqn++)
   {
     // allocate A coefficients for conduction operator
     m_aco[iEqn] =  RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, cellFact));

     // allocate equation diffusion coefficient on regular and irregular cells
     m_bco[iEqn]     = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, fluxFact));
     m_bcoIrreg[iEqn]= RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblgPtr->getDBL(), 1, 4*IntVect::Unit, bivfFact));
   }

  EBAMRRANS* coarPtr = getCoarserLevel();
  EBAMRRANS* finePtr = getFinerLevel();

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

      // maintain levelRANS
      m_ebLevelRANS.define(m_eblgPtr->getDBL(),
                           coEBLGPtr->getDBL(),
                           m_eblgPtr->getEBISL(),
                           coEBLGPtr->getEBISL(),
                           m_eblgPtr->getDomain(),
                           nRefCrse,
                           m_dx,
                           m_hasCoarser,
                           m_hasFiner,
                           m_ebPatchModel);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if(!s_noEBCF)
        {
          CH_TIME("fineToCoar_defs");
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

      // skipping the quadCFI stuff. The smarter me has figured out the reason in the past and successfully solved transport eqn. Im just trusting him

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_ebLevelRANS.define(m_eblgPtr->getDBL(),
                           DisjointBoxLayout(),
                           m_eblgPtr->getEBISL(),
                           EBISLayout(),
                           m_eblgPtr->getDomain(),
                           m_ref_ratio,
                           m_dx,
                           m_hasCoarser,
                           m_hasFiner,
                           m_ebPatchModel);
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
/*********/
void EBAMRRANS::
preRegrid(int                         a_base_level,
          const Vector<Vector<Box> >& a_new_grids)
{

}
/*********/
void EBAMRRANS::
postRegrid(int a_base_level)
{
  if(m_level == a_base_level) defineSolvers();
}
/**********/
void EBAMRRANS::
defineSolvers()
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
     Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);

      EBAMRRANS* coarsestLevel = dynamic_cast<EBAMRRANS*>(hierarchy[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblgPtr->getDomain();
      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRRANS* RANSLevel = (EBAMRRANS*)(hierarchy[ilev]);
          eblgs       [ilev] = *RANSLevel->m_eblgPtr;
          grids       [ilev] = RANSLevel->m_eblgPtr->getDBL();
          refRat      [ilev] = RANSLevel->m_ref_ratio;
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

      Real bottomCushion = 1.0;
      if(pp.contains("bottom_cushion"))
        {
          pp.get("bottom_cushion", bottomCushion);
        }

      s_diffuseSolver.resize(m_nEqn);
      s_diffuseIntegratorBE.resize(m_nEqn);
      s_diffuseIntegratorTGA.resize(m_nEqn);

      for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
       {
         s_diffuseSolver[iEqn] = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());

         s_diffuseSolver[iEqn]->define(lev0Dom, *s_diffuseOpFactory[iEqn], &s_botSolver, nlevels);

         s_diffuseSolver[iEqn]->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);

         s_diffuseSolver[iEqn]->m_verbosity = mgverb; 
         s_diffuseSolver[iEqn]->m_bottomSolverEpsCushion = bottomCushion;

         s_diffuseIntegratorBE[iEqn]  = RefCountedPtr< EBLevelBackwardEuler>( new  EBLevelBackwardEuler(grids, refRat, lev0Dom, s_diffuseOpFactory[iEqn], s_diffuseSolver[iEqn]));
         s_diffuseIntegratorBE[iEqn]->setEBLG(eblgs); 
       }
    }
}
/*********/
void EBAMRRANS:::defineFactories(bool a_atHalfTime)
{
  if (m_params.m_doDiffusion)
    {
      Vector<EBAMRLevel*> hierarchy = EBAMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                             refRat(nlevels);
      Vector<EBLevelGrid>                                     eblgs(nlevels);
      Vector<DisjointBoxLayout>                               grids(nlevels);

      Vector<Vector<RefCountedPtr<LevelData<EBCellFAB> > > >  aco(m_nEqn);
      Vector<Vector<RefCountedPtr<LevelData<EBFluxFAB> > > >  bco(m_nEqn);
      Vector<Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > >    bcoIrreg(m_nEqn);

      Vector<RefCountedPtr<EBQuadCFInterp> >                quadCFI(nlevels);

      for (int i = iEqn; iEqn < m_nEqn; iEqn++)
      {
        aco[iEqn].resize(nlevels);
        bco[iEqn].resize(nlevels);
        bcoIrreg[iEqn].resize(nlevels);
      }

      EBAMRRANS* coarsestLevel = (EBAMRRANS*)(hierarchy[0]);
      Real           lev0Dx      = (coarsestLevel->m_dx[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblgPtr->getDomain();

      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRRANS* RANSLevel = dynamic_cast<EBAMRRANS*>(hierarchy[ilev]);

          EBCellFactory fact(RANSLevel->m_eblgPtr->getEBISL());
          LevelData<EBCellFAB> halfSt(RANSLevel->m_eblgPtr->getDBL(),m_nComp, 4*IntVect::Unit, fact);
          if(a_atHalfTime)
            {
              RANSLevel->getHalfState(halfSt);
            }
          else
            {
              Interval interv(0, m_nComp-1);
              RANSLevel->m_stateOld.copyTo(interv, halfSt, interv);
            }
          RANSLevel->fillCoefficients(halfSt);

          eblgs       [ilev] = *RANSLevel->m_eblgPtr;
          grids       [ilev] = RANSLevel->m_eblgPtr->getDBL();
          refRat      [ilev] = RANSLevel->m_ref_ratio;
          quadCFI     [ilev] = *RANSLevel->m_quadCFIPtr;
 
          for (int iEqn = 0; iEqn < m_nSpec; iEqn++)
           {
             aco[iEqn][ilev] = RANSLevel->m_aco[iEqn];
             bco[iEqn][ilev] = RANSLevel->m_bco[iEqn];
             bcoIrreg[iEqn][ilev] = RANSLevel->m_bcoIrreg[iEqn];  
           } // end eqn loop
        } // end level loop       

      //alpha, beta get replaced in tga solves
      Real alpha = 1;
      Real beta  = 1;
      IntVect  giv    = 4*IntVect::Unit;

      // Diffusion operator.
      int relaxType = 1;

// set the domain and the EB diffusion BCs
// for now 0 neuman
      setBCs();   

      for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
        {
          s_diffuseOpFactory[iEqn] = 
            RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
            (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
            (new EBConductivityOpFactory(eblgs, quadCFI, alpha, beta, aco[iEqn],
                                         bco[iEqn], bcoIrreg[iEqn], lev0Dx,
                                         refRat, m_domBC[iEqn],
                                         m_ebBC[iEqn], giv, giv, relaxType)));  
        }
    }

  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
    {
      m_bco[iEqn]->exchange(Interval(0,0));
    }
}
/*********/
void EBAMRRANS::getHalfState(LevelData<EBCellFAB>& a_stateInt)
{
  //interpolate state to n+1/2
  Real told = 0; Real tnew = 1; Real time = 0.5;
  EBArith::timeInterpolate(a_stateInt, m_stateOld, m_stateNew,
                           m_eblgPtr->getDBL(), time, told, tnew);
}
/*********/
void EBAMRRANS::
fillCoefficients(const LevelData<EBCellFAB>& a_state)
{
// the ugly stuff  
  EBCellFactory cellFact(m_eblgPtr->getEBISL());
  EBFluxFactory fluxFact(m_eblgPtr->getEBISL());
  int nghost = 4;

  LevelData<EBCellFAB> stateCell(m_eblgPtr->getDBL(), m_nComp, nghost*IntVect::Unit, cellFact);
  LevelData<EBFluxFAB> stateFace(m_eblgPtr->getDBL(), m_nComp, IntVect::Zero, fluxFact);

  a_state.copyTo(Interval(0,m_nComp-1), stateCell, Interval(0,m_nComp-1));

  if (m_hasCoarser)
   {
     EBAMRRANS* coarPtr = getCoarserLevel();
     EBPWLFillPatch patcher(m_eblgPtr->getDBL(), coarPtr->m_eblgPtr->getDBL(),
                             m_eblgPtr->getEBISL(), coarPtr->m_eblgPtr->getEBISL(),
                             coarPtr->m_eblgPtr->getDomain(), m_ref_ratio, m_nComp, nghost);

     EBCellFactory coarFact(coarPtr->m_eblgPtr->getEBISL());
     LevelData<EBCellFAB> coarState(coarPtr->m_eblgPtr->getDBL(), m_ncons, nghost*IntVect::Unit, coarFact);

     coarPtr->m_stateNew.copyTo(Interval(0,m_nComp-1), coarState, Interval(0,m_nComp-1));

     patcher.interpolate(stateCell, coarState, coarState, m_time, m_time, m_time, Interval(0,m_nCons-1)); 
   }

  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
   {
     EBLevelDataOps::averageCellToFaces(stateFace, stateCell, m_eblgPtr->getDBL(), m_eblgPtr->getEBISL(), m_eblgPtr->getDomain(), iEqn); 
   }

  setDiffusionCoefficients(stateCell, stateFace);
}
/*********/
void EBAMRRANS::setDiffusionCoefficients(LevelData<EBCellFAB>& a_stateCell,
                                         LevelData<EBFluxFAB>& a_stateFace)
{
  EBFluxFactory fluxFact(m_eblgPtr->getEBISL());
  BaseIVFAB<Real> bivFact(m_eblgPtr->getEBISL(), m_sets);

  LevelData<EBFluxFAB> diffCoeff(m_eblgPtr->getDBL(), m_nEqn, 4*IntVect::Unit, fluxFact);
  LevelData<BaseIVFAB<Real> > diffCoeffIrreg(m_eblgPtr->getDBL(), m_nComp, 4*IntVect::Unit, bivfFact);

  for (int iEqn = 0; iEqn < m_nEqn ; iEqn++)
   {
     m_bco[iEqn]->copyTo(Interval(0,0), diffCoeff, Interval(iEqn,iEqn));
     m_bcoIrreg[iEqn]->copyTo(Interval(0,0), diffCoeffIrreg, Interval(iEqn,iEqn));
   }

  m_ebLevelRANS.getDiffusionCoefficients(diffCoeff, diffCoeffIrreg, a_stateCell, a_stateFace, m_time, m_dt);

  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
   {
     diffCoeff.copyTo(Interval(iEqn,iEqn), *m_bco[iEqn], Interval(0,0));
     diffCoeffIrreg.copyTo(Interval(iEqn,iEqn), *m_bcoIrreg[iEqn], Interval(0,0));
   }
}/********/
void EBAMRRANS::
setBCs()
{
  m_domBC.resize(m_nEqn);
  m_ebBC.resize(m_nEqn);

  // no diffusion at the domain and EB for now
  NeumannConductivityEBBCFactory* ebBCFact =  new NeumannConductivityEBBCFactory();
  ebBCFact->setValue(0.);

  NeumannConductivityDomainBCFactory* domBCFact = new NeumannConductivityDomainBCFactory();
  domBCFact->setValue(0.);

  for (int iEqn = 0; iEqn < m_nEqn; iEqn++)
   {
     m_domBC[iEqn] = RefCountedPtr<BaseDomainBCFactory>(domBCFact);
     m_ebBC[iEqn] = RefCountedPtr<BaseEBBCFactory>(ebBCFact);
   }
}
/*********/
#endif
#include "NamespaceFooter.H"
