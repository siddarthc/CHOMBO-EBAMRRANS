#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "InflowOutflowIBC.H"
#include "InflowOutflowAdvectBC.H"
#include "InflowOutflowPoissonDomainBC.H"
#include "VoFIterator.H"
#include "NeumannPoissonEBBC.H"
#include "InflowOutflowIBC.H"
#include "EBLevelDataOps.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

/***/
RefCountedPtr<BaseEBBCFactory>
InflowOutflowIBC::getPressureEBBC()  const
{
  NeumannPoissonEBBCFactory* neumptr = new NeumannPoissonEBBCFactory();

  neumptr->setValue(0.);

  RefCountedPtr<BaseEBBCFactory>  retval = RefCountedPtr<BaseEBBCFactory>(neumptr);
  return retval;
}
/***/
RefCountedPtr<BaseDomainBCFactory>
InflowOutflowIBC::getVelBC(int a_icomp) const
{
  //this is only used in the helmholtz solve so always no slip
  RefCountedPtr<BaseDomainBCFactory> helmBC = RefCountedPtr<BaseDomainBCFactory>(new InflowOutflowHelmholtzDomainBCFactory(m_flowDir,

                                           m_doJet2,

                                           m_jet1inflowVel,

                                           m_jet2inflowVel,

                                           m_doJet1PoiseInflow,

                                           m_doJet2PoiseInflow,

                                           m_doSlipWallsHi,

                                           m_doSlipWallsLo,

                                           m_jet1PoiseInflowFunc,

                                           m_jet2PoiseInflowFunc));  
  return helmBC;
}

///
RefCountedPtr<BaseEBBCFactory>
InflowOutflowIBC::getVelocityEBBC(int a_velComp) const
{
  RefCountedPtr<BaseEBBCFactory>  retval;
  //slip walls only apply to the domain boundary
  DirichletPoissonEBBCFactory* dirptr = new DirichletPoissonEBBCFactory();
  Real val = 0.0;
  dirptr->setValue(val);
  dirptr->setOrder(m_orderEBBC);
  retval = RefCountedPtr<BaseEBBCFactory>(dirptr);

  return retval;
}
///
RefCountedPtr<BaseDomainBCFactory>
InflowOutflowIBC::
getPressBC() const
{
  RefCountedPtr<BaseDomainBCFactory> poisBC
    = RefCountedPtr<BaseDomainBCFactory>
        (new InflowOutflowPoissonDomainBCFactory(m_flowDir,
                                                 m_doJet2,
                                                 m_jet1inflowVel,
                                                 m_jet2inflowVel,
                                                 m_doJet1PoiseInflow,
                                                 m_doJet2PoiseInflow,
                                                 m_doSlipWallsHi,
                                                 m_doSlipWallsLo,
                                                 m_jet1PoiseInflowFunc,
                                                 m_jet2PoiseInflowFunc));
  return poisBC;
}
///
RefCountedPtr<EBPhysIBCFactory>
InflowOutflowIBC::
getVelAdvectBC(int a_velComp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval
    = RefCountedPtr<EBPhysIBCFactory>
        (new InflowOutflowAdvectBCFactory(m_flowDir,
                                          m_doJet2,
                                          m_jet1inflowVel,
                                          m_jet2inflowVel,
                                          a_velComp,
                                          m_doJet1PoiseInflow,
                                          m_doJet2PoiseInflow,
                                          m_doSlipWallsHi,
                                          m_doSlipWallsLo,
                                          m_jet1PoiseInflowFunc,
                                          m_jet2PoiseInflowFunc));
  return retval;
}
///
RefCountedPtr<EBPhysIBCFactory>
InflowOutflowIBC::
getScalarAdvectBC(const int&  a_comp) const
{
  RefCountedPtr<EBPhysIBCFactory> retval = RefCountedPtr<EBPhysIBCFactory>(new ExtrapAdvectBCFactory());
  return retval;
}

///
RefCountedPtr<BaseDomainBCFactory>
InflowOutflowIBC::
getMACVelBC() const
{
  RefCountedPtr<BaseDomainBCFactory> poisBC
    = RefCountedPtr<BaseDomainBCFactory>
        (new InflowOutflowPoissonDomainBCFactory(m_flowDir,
                                                 m_doJet2,
                                                 m_jet1inflowVel,
                                                 m_jet2inflowVel,
                                                 m_doJet1PoiseInflow,
                                                 m_doJet2PoiseInflow,
                                                 m_doSlipWallsHi,
                                                 m_doSlipWallsLo,
                                                 m_jet1PoiseInflowFunc,
                                                 m_jet2PoiseInflowFunc));
  return poisBC;
}

///
void
InflowOutflowIBC::
initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
   EBLevelDataOps::setVal(a_velocity, 0.0, m_flowDir);
  if (m_initPoiseData)
    {
//      MayDay::Error("initPoiseData not set for Velocity");
      CH_assert(m_isFortranCommonSet);
      CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);
      Real dx = a_dx[0];
      // Iterator of all grids in this level
      for (DataIterator dit = a_velocity.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = a_ebisl[dit()];
          // Storage for current grid
          EBCellFAB& velFAB = a_velocity[dit()];
    
          BaseFab<Real>& regVelFAB = velFAB.getSingleValuedFAB();
          // Box of current grid
          Box uBox = regVelFAB.box();
          uBox &= a_domain;
          // Set up initial condition in this grid
          /**/
          FORT_JET2VELINIT(CHF_CONST_FRA(regVelFAB),
                           CHF_CONST_REAL(dx),
                           CHF_BOX(uBox));
          /**/

          //now for the multivalued cells.  Since it is pointwise,
          //the regular calc is correct for all single-valued cells.
          IntVectSet ivs = ebisBox.getMultiCells(uBox);
          for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect& iv = vof.gridIndex();
              RealVect velocity;
              /**/
              FORT_POINTJET2VELINIT(CHF_REALVECT(velocity),
                                    CHF_CONST_INTVECT(iv),
                                    CHF_CONST_REAL(dx));
              /**/    
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  velFAB(vof, idir) = velocity[idir];
//                  velFAB(vof, idir) = 0.; 
                }

            }//end loop over multivalued cells
       } //end loop over boxes 
    }
}
///
void
InflowOutflowIBC::
initializePressureGradient(LevelData<EBCellFAB>&    a_gradient,
                           const DisjointBoxLayout& a_grids,
                           const EBISLayout&        a_ebisl,
                           const ProblemDomain&     a_domain,
                           const RealVect&          a_origin,
                           const Real&              a_time,
                           const RealVect&          a_dx) const
{
  EBLevelDataOps::setToZero(a_gradient);
  if (m_initPoiseData)
    {
//      MayDay::Error("initPoiseData not set for pressure");

      CH_assert(m_isFortranCommonSet);
      CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);
      Real dx = a_dx[0];
      // Iterator of all grids in this level
      for (DataIterator dit = a_gradient.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = a_ebisl[dit()];
          // Storage for current grid
          EBCellFAB& gradFAB = a_gradient[dit()];
    
          BaseFab<Real>& regGradFAB = gradFAB.getSingleValuedFAB();
          // Box of current grid
          Box uBox = regGradFAB.box();
          uBox &= a_domain;
          // Set up initial condition in this grid
          /**/
          FORT_JET2GRADPINIT(CHF_CONST_FRA(regGradFAB),
                             CHF_CONST_REAL(dx),
                             CHF_BOX(uBox));
          /**/

          //now for the multivalued cells.  Since it is pointwise,
          //the regular calc is correct for all single-valued cells.
          IntVectSet ivs = ebisBox.getMultiCells(uBox);
          for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect& iv = vof.gridIndex();
              RealVect grad;
              /**/
              FORT_POINTJET2GRADPINIT(CHF_REALVECT(grad),
                                      CHF_CONST_INTVECT(iv),
                                      CHF_CONST_REAL(dx));
              /**/    
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  gradFAB(vof, idir) = grad[idir];
                }

            }//end loop over multivalued cells

       } //end loop over boxes 
    }
}
///
void
InflowOutflowIBC::
initializePressure(LevelData<EBCellFAB>&    a_pressure,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  EBLevelDataOps::setToZero(a_pressure);
}
///
void
InflowOutflowIBC::
initializeScalar ( LevelData<EBCellFAB>&    a_scalar,
                   const DisjointBoxLayout& a_grids,
                   const EBISLayout&        a_ebisl,
                   const ProblemDomain&     a_domain,
                   const RealVect&          a_origin,
                   const Real&              a_time,
                   const RealVect&          a_dx) const
{
  EBLevelDataOps::setToZero(a_scalar);
}

#include "NamespaceFooter.H"
