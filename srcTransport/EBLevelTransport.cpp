#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DebugOut.H"
#include "EBLevelTransport.H"
#include "BaseIVFactory.H"
#include "BaseIFFactory.H"
#include "BaseIFFAB.H"
#include "EBFluxFAB.H"
#include "FaceIterator.H"
#include "REAL.H"
#include "EBCellFactory.H"
#include "FaceIterator.H"
#include "EBPatchTransportF_F.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "EBAMRIO.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

/*****************************/
EBLevelTransport::
EBLevelTransport()
{
  m_isDefined = false;
}
/*****************************/
bool
EBLevelTransport::
isDefined() const
{
  return m_isDefined;
}
/*****************************/
EBLevelTransport::
~EBLevelTransport()
{
}
/*****************************/
EBLevelTransport::
EBLevelTransport(const DisjointBoxLayout&           a_thisDBL,
                 const DisjointBoxLayout&           a_coarDBL,
                 const EBISLayout&                  a_thisEBISL,
                 const EBISLayout&                  a_coarEBISL,
                 const ProblemDomain&               a_domain,
                 const int&                         a_nRefine,
                 const RealVect&                    a_dx,
                 const bool&                        a_hasCoarser,
                 const bool&                        a_hasFiner,
                 RefCountedPtr<EBPatchTransport>&   a_patchGodunov,
                 const EBIndexSpace*          const a_eb)
{
  CH_TIME("EBLevelTransport::EBLevelTransport");
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
         a_patchGodunov,
         a_eb);
}

/*****************************/
void
EBLevelTransport::
define(const DisjointBoxLayout&           a_thisDBL,
       const DisjointBoxLayout&           a_coarDBL,
       const EBISLayout&                  a_thisEBISL,
       const EBISLayout&                  a_coarEBISL,
       const ProblemDomain&               a_domain,
       const int&                         a_nRefine,
       const RealVect&                    a_dx,
       const bool&                        a_hasCoarser,
       const bool&                        a_hasFiner,
       RefCountedPtr<EBPatchTransport>&     a_patchGodunov,
       const EBIndexSpace*          const a_eb)
{
  CH_TIME("EBLevelTransport::define");
  CH_assert(a_dx[0] > 0.0);
  CH_assert(a_nRefine > 0);

  m_ebPatchGodunov = a_patchGodunov;

  m_isDefined = true;
  m_thisGrids    = a_thisDBL;
  m_thisEBISL    = a_thisEBISL;
  m_refRatCrse   = a_nRefine;
  m_dx           = a_dx;
  m_domain       = a_domain;
  m_hasCoarser   = a_hasCoarser;
  m_hasFiner     = a_hasFiner;
  if (m_hasCoarser)
    {
      m_coarGrids = a_coarDBL;
      m_coarEBISL = a_coarEBISL;
    }

  m_nVar         =  1;
  m_nPrim        =  1;
  m_nCons        =  1;
  m_nFlux        =  1;

  m_irregSetsSmall.define(m_thisGrids);
  for(DataIterator dit = m_thisGrids.dataIterator();
      dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if(!ebisBox.isAllCovered())
        {
          const Box&     thisBox = m_thisGrids.get(dit());
          m_irregSetsSmall[dit()] = ebisBox.getIrregIVS(thisBox);

        }
    }

  m_cfIVS.define(m_thisGrids); 
  for(DataIterator ditOuter = m_thisGrids.dataIterator();
      ditOuter.ok(); ++ditOuter)
    {
      const EBISBox& ebisBox = m_thisEBISL[ditOuter()];
      if(!ebisBox.isAllCovered())
        {
          // make a grown box of the grid and then subtract off each grid in the
          //domain to make the coarse-fine IVS
          // CFIVS = grow(b, 1) - sum(bi)
          Box grownBox = grow(m_thisGrids.get(ditOuter()), 1);
          grownBox &= m_domain;
          IntVectSet complementIVS(grownBox);

          for(LayoutIterator litInner = m_thisGrids.layoutIterator();
              litInner.ok(); ++litInner)
            {
              Box innerBox = m_thisGrids.get(litInner());
              complementIVS -= innerBox;
            }

          m_cfIVS[ditOuter()] = complementIVS;
        }
    }

  m_nGhost = a_thisEBISL.getGhost();

  if (m_hasCoarser)
    {
      CH_TIME("EBLevelTransport::define::fillPatchDefine");
      ProblemDomain domainCrse = coarsen(m_domain, m_refRatCrse);

      //patcher is defined with the number of conserved vars.
      m_fillPatch.define(m_thisGrids, m_coarGrids,
                         m_thisEBISL, m_coarEBISL,
                         domainCrse, m_refRatCrse, m_nVar,
                         m_nGhost, a_eb);

      m_fillPatchVel.define(m_thisGrids, m_coarGrids,
                            m_thisEBISL, m_coarEBISL,
                            domainCrse, m_refRatCrse, SpaceDim,
                            m_nGhost, a_eb);

    }

  m_irregSetsSmall.define(m_thisGrids);
  for(DataIterator dit = m_thisGrids.dataIterator();
      dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if(!ebisBox.isAllCovered())
        {
          const Box&     thisBox = m_thisGrids.get(dit());
          m_irregSetsSmall[dit()] = ebisBox.getIrregIVS(thisBox);

        }
    }
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      EBArith::defineFluxInterpolant(m_fluxInterpolants[faceDir],
                                     m_irregSetsGrown  [faceDir],
                                     m_thisGrids, m_thisEBISL,
                                     m_domain, m_nFlux, faceDir);
    }

  BaseIVFactory<Real> cellFactorySmall(m_thisEBISL, m_irregSetsSmall);
  m_nonConsDivergence.define(m_thisGrids, m_nCons,
                             IntVect::Zero, cellFactorySmall);
  m_ebIrregFaceFlux.define(m_thisGrids, m_nCons,
                           IntVect::Zero, cellFactorySmall);

 //create temp data with the correct number of ghost cells
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  EBCellFactory factory(m_thisEBISL);
  IntVect flatGhostIV = 3*IntVect::Unit;
  m_flattening.define(m_thisGrids, 1, flatGhostIV, factory); 
}
/*****************************/
void EBLevelTransport::
floorPrimitives(LevelData<EBCellFAB>&         a_primState, 
                Real a_time, Real a_dt)
{
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& primState = a_primState[dit()];
      const IntVectSet& cfivs = m_cfIVS[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      const Box& cellBox = m_thisGrids.get(dit());
      m_ebPatchGodunov->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
      m_ebPatchGodunov->floorPrimitives(primState, cellBox);
    }
}
/*****************************/
void EBLevelTransport::
floorConserved(LevelData<EBCellFAB>&         a_consState,
                Real a_time, Real a_dt)
{ 
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& consState = a_consState[dit()];
      const IntVectSet& cfivs = m_cfIVS[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      const Box& cellBox = m_thisGrids.get(dit());
      m_ebPatchGodunov->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
      m_ebPatchGodunov->floorConserved(consState, cellBox);
    }
}
/*****************************/
void EBLevelTransport::
divergeF(LevelData<EBCellFAB>&         a_divergeF,
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
         const Real&                   a_dt)
{
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  CH_assert(isDefined());

  CH_assert(a_consState.disjointBoxLayout() == m_thisGrids);

  //fill a_consState with data interpolated between constate and coarser data
  fillCons(a_consState, a_normalVel, a_consStateCoarseOld, a_consStateCoarseNew, a_normalVelCoarseOld, a_normalVelCoarseNew, a_time, a_coarTimeOld, a_coarTimeNew); 

  // clear flux registers with fine level
  //remember this level is the coarse level to the fine FR
  if(m_hasFiner)
    {
      a_fineFluxRegister.setToZero();
    }

  //compute flattening coefficients. this saves a ghost cell. woo hoo.
  computeFlattening(a_consState, a_time, a_dt);

  //this includes copying flux into flux interpolant and updating
  //regular grids and incrementing flux registers.

  doRegularUpdate(a_divergeF, a_consState, a_normalVel, a_advVel, a_coveredAdvVelMinu, a_coveredAdvVelPlus, a_fineFluxRegister, a_coarFluxRegister, a_source, a_time, a_dt);

  //this does irregular update and deals with flux registers.
  //also computes the mass increment
  doIrregularUpdate(a_divergeF, a_consState, a_fineFluxRegister, a_coarFluxRegister, a_massDiff, a_time, a_dt);

}
/*****************************/
void
EBLevelTransport::
fillCons(LevelData<EBCellFAB>&         a_consState,
         LevelData<EBCellFAB>&         a_normalVel,
         const LevelData<EBCellFAB>&   a_consStateCoarseOld,
         const LevelData<EBCellFAB>&   a_consStateCoarseNew,
         const LevelData<EBCellFAB>&   a_normalVelCoarseOld,
         const LevelData<EBCellFAB>&   a_normalVelCoarseNew,
         const Real&                   a_time,
         const Real&                   a_coarTimeOld,
         const Real&                   a_coarTimeNew)
{
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  Interval intervSD(0,SpaceDim-1);

  if (m_hasCoarser)
    {
      m_fillPatch.interpolate(a_consState,
                              a_consStateCoarseOld,
                              a_consStateCoarseNew,
                              a_coarTimeOld,
                              a_coarTimeNew,
                              a_time,
                              consInterv);

      m_fillPatchVel.interpolate(a_normalVel,
                                 a_normalVelCoarseOld,
                                 a_normalVelCoarseNew,
                                 a_coarTimeOld,
                                 a_coarTimeNew,
                                 a_time,
                                 intervSD);
    }

  // Exchange all the data between grids
  a_consState.exchange(consInterv);
  a_normalVel.exchange(intervSD);
}
/***************/
void
EBLevelTransport::
computeFlattening(LevelData<EBCellFAB>& a_cons, Real a_time, Real a_dt)
{
  //compute flattening coefficients. this saves a ghost cell. woo hoo.
  int ibox = 0;
  bool verbose = false;
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& consState = a_cons[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if(!ebisBox.isAllCovered())
        {
          const Box& cellBox = m_thisGrids.get(dit());
          const IntVectSet& cfivs = m_cfIVS[dit()];
          m_ebPatchGodunov->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
          m_ebPatchGodunov->setCoveredConsVals(consState);
          int nPrim  = m_ebPatchGodunov->numPrimitives();
          EBCellFAB primState(ebisBox, consState.getRegion(), nPrim);
          int logflag = 0;
          //debug
          //verbose = true;
          //end debug
          m_ebPatchGodunov->consToPrim(primState, consState, consState.getRegion(), logflag, verbose);
          EBCellFAB& flatteningFAB = m_flattening[dit()];
          //this will set the stuff over the coarse-fine interface
          flatteningFAB.setVal(1.);
          if(m_ebPatchGodunov->usesFlattening())
            {
              m_ebPatchGodunov->computeFlattening(flatteningFAB,
                                                  primState,
                                                  cellBox);
            }
          ibox++;
        }
    }
  Interval zerointerv(0,0);
  m_flattening.exchange(zerointerv);
}
/***************/
void
EBLevelTransport::
doRegularUpdate(LevelData<EBCellFAB>&         a_divF,
                LevelData<EBCellFAB>&         a_cons,
                const LevelData<EBCellFAB>&   a_normalVel,
                const LevelData<EBFluxFAB>&   a_advVel,
                const LayoutData< Vector <BaseIVFAB<Real> * > >& a_coveredAdvVelMinu,
                const LayoutData< Vector <BaseIVFAB<Real> * > >& a_coveredAdvVelPlus,
                EBFluxRegister&               a_fineFluxRegister,
                EBFluxRegister&               a_coarFluxRegister,
                const LevelData<EBCellFAB>&   a_source,
                Real a_time, Real a_dt)
{
  bool verbose = false;
  int ibox = 0;
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);

    for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      const Box& cellBox = m_thisGrids.get(dit());
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if(!ebisBox.isAllCovered())
        {
          const IntVectSet& cfivs = m_cfIVS[dit()];

          EBCellFAB& consState = a_cons[dit()];
          m_ebPatchGodunov->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt); 
          //begin debug
          //int here = 0;
          //if(cellBox.contains(EBPatchGodunov::s_debugIV) && (EBPatchGodunov::s_whichLev==1))
          //  {
          //    here = 1;
          //  }
          //end debug

          const EBCellFAB& source = a_source[dit()];
          const EBCellFAB& normalVel = a_normalVel[dit()];
          const EBFluxFAB& advVel = a_advVel[dit()];
          const Vector<BaseIVFAB <Real>* >& coveredAdvVelMinu = a_coveredAdvVelMinu[dit()];
          const Vector<BaseIVFAB <Real>* >& coveredAdvVelPlus = a_coveredAdvVelPlus[dit()];

          m_ebPatchGodunov->setVelocities(normalVel, advVel, coveredAdvVelMinu, coveredAdvVelPlus);

          EBFluxFAB flux(ebisBox, cellBox, m_nFlux);
          BaseIVFAB<Real>& nonConsDiv    = m_nonConsDivergence[dit()];
          BaseIVFAB<Real>& ebIrregFlux   = m_ebIrregFaceFlux[dit()];
          flux.setVal(7.89);
          ebIrregFlux.setVal(7.89);
          const IntVectSet& ivsIrreg     = m_irregSetsSmall[dit()];
          const EBCellFAB& flatteningFAB = m_flattening[dit()];

          BaseIVFAB<Real>  coveredPrimMinu[SpaceDim];
          BaseIVFAB<Real>  coveredPrimPlus[SpaceDim];
          Vector<VolIndex> coveredFaceMinu[SpaceDim];
          Vector<VolIndex> coveredFacePlus[SpaceDim];
          EBFluxFAB facePrim;

          EBCellFAB& divFFAB = a_divF[dit()];
//          EBPatchTransport& patchTrans = (EBPatchTransport&) (*m_ebPatchGodunov);

          EBPatchGodunov::setCurComp(0);
          EBPatchGodunov::setDoingVel(0);
          EBPatchGodunov::setDoingAdvVel(0);


          m_ebPatchGodunov->primitivesAndDivergences(divFFAB, consState,
                                             facePrim,
                                             coveredPrimMinu,
                                             coveredPrimPlus,
                                             coveredFaceMinu,
                                             coveredFacePlus,
                                             flux, ebIrregFlux,
                                             nonConsDiv,flatteningFAB,
                                             source, cellBox, ivsIrreg,
                                             dit(),verbose);

          //do fluxregister cha-cha
          /*
            Coarse flux register is flux register with the coarse level.
            Fine flux register is the flux register with the fine level.
            To the finer level FR, this level is the coarse level.
            To the coarser level FR, this level is the fine level.
          */
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real scale = a_dt;

              EBFaceFAB fluxRegFlux;
              if(m_hasFiner)
                {
                  a_fineFluxRegister.incrementCoarseRegular(flux[idir], scale,dit(),
                                                            consInterv, idir);
                }

              if(m_hasCoarser)
                {
                  for(SideIterator sit; sit.ok(); ++sit)
                    {
                      a_coarFluxRegister.incrementFineRegular(flux[idir],scale, dit(),
                                                              consInterv, idir,sit());
                    }
                }
            }

          //copy fluxes into sparse interpolant
          for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              IntVectSet ivsIrregGrown = m_irregSetsGrown[faceDir][dit()];
              ivsIrregGrown &= cellBox;
              FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

              BaseIFFAB<Real>& interpol = m_fluxInterpolants[faceDir][dit()];
              interpol.setVal(7.7777e7);
              EBFaceFAB& fluxDir = flux[faceDir];
              for(FaceIterator faceit(ivsIrregGrown, ebisBox.getEBGraph(),
                                      faceDir, stopCrit);
                  faceit.ok(); ++faceit)
                {
                  for(int ivar = 0; ivar < m_nFlux; ivar++)
                    {
                      interpol(faceit(), ivar) = fluxDir(faceit(), ivar);
                    }
                }

            }
        }
    }

  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      m_fluxInterpolants[faceDir].exchange(fluxInterv);
    } 
}
/***************/
void
EBLevelTransport::
doIrregularUpdate(LevelData<EBCellFAB>&         a_divergeF,
                  LevelData<EBCellFAB>&         a_cons,
                  EBFluxRegister&               a_fineFluxRegister,
                  EBFluxRegister&               a_coarFluxRegister,
                  LevelData<BaseIVFAB<Real> >&  a_massDiff,
                  Real a_time, Real a_dt)
{
  //now do the irregular update
  int ibox = 0;
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      const Box& cellBox = m_thisGrids.get(dit());
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      if(!ebisBox.isAllCovered())
        {
          const IntVectSet& cfivs = m_cfIVS[dit()];

          EBCellFAB& consState = a_cons[dit()];
          BaseIVFAB<Real>& redMass = a_massDiff[dit()];

          m_ebPatchGodunov->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
          //begin debug
          //int here = 0;
          //if(cellBox.contains(EBPatchGodunov::s_debugIV) && (EBPatchGodunov::s_whichLev==1))
          //  {
          //    here = 1;
          //  }
          //end debug

          BaseIFFAB<Real> centroidFlux[SpaceDim];
          const BaseIFFAB<Real>* interpolantGrid[SpaceDim];
          const IntVectSet& ivsIrregSmall = m_irregSetsSmall[dit()];
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& interpol = m_fluxInterpolants[idir][dit()];
              interpolantGrid[idir] = &interpol;
              BaseIFFAB<Real>& fluxDir= centroidFlux[idir];
              fluxDir.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, m_nFlux);
            }
          m_ebPatchGodunov->interpolateFluxToCentroids(centroidFlux,
                                                       interpolantGrid,
                                                       ivsIrregSmall);

          //update the state and interpolate the flux
          const BaseIVFAB<Real>& nonConsDiv = m_nonConsDivergence[dit()];
          const BaseIVFAB<Real>& ebIrregFlux = m_ebIrregFaceFlux[dit()];

          m_ebPatchGodunov->hybridDivergence(a_divergeF[dit()], consState,  redMass,
                                             centroidFlux, ebIrregFlux, nonConsDiv,
                                             cellBox, ivsIrregSmall);


          //do fluxregister mambo
          /*
            Coarse flux register is flux register with the coarse level.
            Fine flux register is the flux register with the fine level.
            To the finer level FR, this level is the coarse level.
            To the coarser level FR, this level is the fine level.
          */
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real scale = a_dt;

              BaseIFFAB<Real> fluxRegFlux;
              if(m_hasFiner)
                {
                  a_fineFluxRegister.incrementCoarseIrregular(centroidFlux[idir],
                                                              scale,dit(),
                                                              consInterv, idir);
                }

              if(m_hasCoarser)
                {
                  for(SideIterator sit; sit.ok(); ++sit)
                    {
                      a_coarFluxRegister.incrementFineIrregular(centroidFlux[idir],
                                                                scale, dit(),
                                                                consInterv, idir,sit());
                    }
                }
            }
        }
    }// end of loop over grids.
}
/*****************************/
void EBLevelTransport::
computeNormalVel(LevelData<EBCellFAB>&                          a_normalVel,
                 const LevelData<EBFluxFAB>&                    a_advectionVel, 
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloLo,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloHi)
{
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& normalVel = a_normalVel[dit()];
      const IntVectSet& cfivs = m_cfIVS[dit()];
      const EBISBox& ebisBox = m_thisEBISL[dit()];
      const Box& cellBox = m_thisGrids.get(dit());
      m_ebPatchGodunov->setValidBox(cellBox, ebisBox, cfivs, 0.0, 0.0);
      m_ebPatchGodunov->averageVelToCC(normalVel, a_advectionVel[dit()], a_coveredVeloLo[dit()], a_coveredVeloHi[dit()], cellBox);
    }  
}
/*****************************/
#include "NamespaceFooter.H"
