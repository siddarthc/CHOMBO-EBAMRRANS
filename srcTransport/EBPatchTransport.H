#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHTRANSPORT_H_
#define _EBPATCHTRANSPORT_H_

#include "EBPatchGodunov.H"
#include "EBPatchGodunovFactory.H"
#include "NamespaceHeader.H"

#define TOL 1.e-12

///
/**
   An EBPatchGodunov for passive advection with an
   input velocity of a particular function value.
*/
class EBPatchTransport : public EBPatchGodunov
{
public:

  ///
  EBPatchTransport();

  ///
  ~EBPatchTransport();

  ///
  /**
     Return number of components for primitive variables.
  */
  virtual int numPrimitives() const
  {
    return 1;
  }

  ///
  /**
     Returns number of components for flux variables.
  */
  virtual int numFluxes() const
  {
    return 1;
  }

  ///
  /**
     Returns number of components for conserved variables.
  */
  virtual int numConserved() const
  {
    return 1;
  }

  ///
  /**
   */
  virtual int numSlopes() const
  {
    return 1;
  }

  virtual void
  slope(EBCellFAB&       a_slopePrim,
        EBCellFAB&       a_slopeNLim,
        const EBCellFAB& a_primState,
        const EBCellFAB& a_flattening,
        const int&       a_dir,
        const Box&       a_box,
        bool a_doAgg = false); 

  virtual void
  doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                      EBCellFAB&       a_deltaWL,
                      EBCellFAB&       a_deltaWR,
                      EBCellFAB&       a_deltaWC,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const Box&       a_box);

  ///
  /**
   */
  void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box);

  ///
  /**
   */

  /// floors if m_isMaxMinSet
  virtual void
  floorPrimitives(EBCellFAB&  a_primState,
                  const Box&  a_box);

  /// floors if m_isMaxMinSet
  virtual void
  floorPrimitives(BaseIVFAB<Real>&   a_primState,
                  const IntVectSet&  a_ivsIrreg);

  /// floors if m_isMaxMinSet
  virtual void
  floorConserved(EBCellFAB&  a_consState,
                 const Box&  a_box);

  /// floors if m_isMaxMinSet
  virtual void
  floorConserved(BaseIVFAB<Real>&   a_consState,
                 const IntVectSet&  a_ivsIrreg);

  /// default noop
  virtual void  getCoveredValuesPrim(Vector<Real>& a_covValues)
  {;}

  ///default noop
  virtual void  getCoveredValuesCons(Vector<Real>& a_covValues)
  {;}

  ///default noop
  virtual void setCoveredConsVals(EBCellFAB& a_consState);

  ///
  /**
     just a copy as default
  */
  virtual void
  consToPrim(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box,
             int              a_logflag,
             bool             a_verbose=false);

  ///
  /**
     just a copy as default
   */
  virtual void
  consToPrim(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

  ///
  /**
     just a copy
   */
  virtual void
  consToPrim(BaseIVFAB<Real>&  a_primState,
             const EBCellFAB&  a_conState,
             const IntVectSet& a_ivs);

  ///
  /**
     just a copy
   */
  virtual void
  primToCons(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box);

  ///
  /**
     just a copy
   */
  virtual void
  primToCons(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

  virtual void extrapolatePrim(EBFluxFAB&        a_facePrim,
                EBCellFAB&                       a_primState,
                EBCellFAB                        a_slopePrim[SpaceDim],
                EBCellFAB                        a_slopeNLim[SpaceDim],
                BaseIVFAB<Real>                  a_coveredPrimMinu[SpaceDim],
                BaseIVFAB<Real>                  a_coveredPrimPlus[SpaceDim],
                Vector<VolIndex>                 a_coveredFaceMinu[SpaceDim],
                Vector<VolIndex>                 a_coveredFacePlus[SpaceDim],
                const EBCellFAB&                 a_flattening,
                const EBCellFAB&                 a_consState,
                const EBCellFAB&                 a_source,
                const Box&                       a_box,
                const DataIndex&                 a_dit,
                bool                             a_verbose);

  ///
  /**
     sends back primitive state (NOT FLUX) using upwinding
  */
  virtual void
  riemann(EBFaceFAB&       a_primGdnv,
          const EBCellFAB& a_primLeft,
          const EBCellFAB& a_primRight,
          const int&       a_dir,
          const Box&       a_box);

  ///
  /**
     sends back primitive state (NOT FLUX) using upwinding
  */
  virtual void
  riemann(BaseIVFAB<Real>&        a_coveredPrim,
          const BaseIVFAB<Real>&  a_extendedState,
          const EBCellFAB&        a_primState,
          const Vector<VolIndex>& a_region,
          const int&              a_dir,
          const Side::LoHiSide&   a_sd,
          const Box& a_box);

  ///
  /**
     sets eb flux to zero
   */
  virtual void
  computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_source);

  ///
  /**
   */
  //Things you do not want to mess with unless you really know what you are doing
  ///
  /**
     always false for now
  */
  virtual bool usesFlattening() const;

  ///
  /**
     always false for now
  */
  virtual bool usesArtificialViscosity() const;


  // virtual void
  // coveredExtrapSlopes(Real&            a_dqc,
  //                     const VolIndex&  a_vof,
  //                     const EBCellFAB& a_primState,
  //                     const int&       a_dir,
  //                     const int&       a_ivar);

  // virtual void
  // coveredExtrapSlopesEBPA(Real&            a_dqc,
  //                         const VolIndex&  a_vof,
  //                         const EBCellFAB& a_primState,
  //                         const Real&      a_wc,
  //                         const int&       a_dir,
  //                         const int&       a_ivar);


  ///
  /**
     Return true if you are using fourth-order slopes.
     Return false if you are using second-order slopes.
  */
  virtual bool usesFourthOrderSlopes() const;

  ///
  /**
     default does -= u_idir d cons/ dx_idir.
     assumes fluxes are really prims.
  */
  virtual void
  updatePrim(EBCellFAB&              a_primMinu,
             EBCellFAB&              a_primPlus,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale); 

  virtual void setVelocities(const EBCellFAB& a_normalVel,
                             const EBFluxFAB& a_advVel,
                             const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelMinu,
                             const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelPlus);

  virtual void setAdvVel(const EBFluxFAB& a_advVel,
                         const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelMinu,
                         const Vector <BaseIVFAB <Real>* >& a_coveredAdvVelPlus);

  virtual void setNormalVel(const EBCellFAB& a_normalVel);

  void
  averageVelToCC(EBCellFAB&                        a_normalVel,
                 const EBFluxFAB&                  a_advectionVel,
                 const Vector<BaseIVFAB<Real> * >& a_coveredVeloMinu,
                 const Vector<BaseIVFAB<Real> * >& a_coveredVeloPlus,
                 const Box&                        a_box) const;

  virtual void
  setValidBox(const Box& a_validBox,
              const EBISBox& a_ebisBox,
              const IntVectSet& a_coarseFineIVS,
              const Real& a_time,
              const Real& a_dt);

  virtual void setTimeAndDt(const Real& a_time, const Real& a_dt)
  {
    m_time = a_time;
    m_dt   = a_dt;
  }

  ///sets generic defaults
  virtual Vector<string> stateNames();
  ///sets generic defaults
  virtual Vector<string> primNames();

  ///
  void useLimiting(bool a_limiting);

  // Things that should not be called because they do not apply

  /// should not be called
  Real
  getMaxWaveSpeed(const EBCellFAB& a_consState,
                  const Box& a_box)
  {
    MayDay::Error("should not be called ");
    return -1.0;
  }
  /// should not be called
  void
  getFlux(EBFluxFAB&       a_flux,
          const EBFluxFAB& a_prim);
  ///
  /**
     should not be called
  */
  Interval velocityInterval() const
  {
    MayDay::Error("should not be called");
    Interval retval;
    return retval;
  }

  ///
  /**
     only used as a comparison with ivar for flattening
     at covered face extrapolation
     don't want flattening-type stuff to pure advected vars
     also used for tagging cells.
  */
  virtual int densityIndex() const
  {
    return 0;
  }


  ///
  /**
     only used as a comparison with ivar for flattening
     at covered face extrapolation
     don't want flattening-type stuff to pure advected vars
  */
  int pressureIndex() const
  {
    return -1;
  }

  ///     should not be called
  int bulkModulusIndex() const
  {
    MayDay::Error("should not be called");
    return -1;
  }

  Real artificialViscosityCoefficient() const;

  void setMaxMin(const Real& a_maxVal,
                 const Real& a_minVal)
  {
    m_isMaxMinSet = true;
    m_maxVal = a_maxVal;
    m_minVal = a_minVal;
  }
  EBCellFAB& getPrimState()
  {
    return m_primState;
  }

  void
  primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                           EBCellFAB&          a_consState,
                           EBFluxFAB&          a_facePrim,
                           BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                           BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                           Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                           Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                           EBFluxFAB&          a_flux,
                           BaseIVFAB<Real>&    a_ebIrregFlux,
                           BaseIVFAB<Real>&    a_nonConservativeDivergence,
                           const EBCellFAB&    a_flattening,
                           const EBCellFAB&    a_source,
                           const Box&          a_box,
                           const IntVectSet&   a_ivsSmall,
                           const DataIndex&    a_dit,
                           bool                a_verbose);

  virtual void
  advectiveDerivative(EBCellFAB&                      a_udotDelRho,
                      const EBFluxFAB&                a_faceRho,
                      const EBFluxFAB&                a_faceVel,
                      const BaseIVFAB<Real>           a_coveredRhoMinu[SpaceDim],
                      const BaseIVFAB<Real>           a_coveredRhoPlus[SpaceDim],
                      const Vector<BaseIVFAB<Real>*>& a_coveredVelMinu,
                      const Vector<BaseIVFAB<Real>*>& a_coveredVelPlus,
                      const Vector<VolIndex>          a_coveredFaceMinu[SpaceDim],
                      const Vector<VolIndex>          a_coveredFacePlus[SpaceDim],
                      const Box&                      a_box);

  bool checkCoveredFaceStencil(const VolIndex& a_vof,
                               const Vector<VolIndex>& a_vofsStencil,
                               const RealVect& a_normal,
                               const int&      a_faceDir,
                               const Side::LoHiSide& a_sd);
  bool checkCoveredFaceStencil2d(const VolIndex& a_vof,
                                 const Vector<VolIndex>& a_vofsStencil,
                                 const RealVect& a_normal,
                                 const int&      a_faceDir,
                                 const Side::LoHiSide& a_sd);
  bool checkCoveredFaceStencil3d(const VolIndex& a_vof,
                                 const Vector<VolIndex>& a_vofsStencil,
                                 const RealVect& a_normal,
                                 const int&      a_faceDir,
                                 const Side::LoHiSide& a_sd);

protected:

  const EBCellFAB* m_normalVelPtr;
  const EBFluxFAB* m_advVelPtr;
  const Vector< BaseIVFAB< Real>* >* m_coveredAdvVelMinuPtr;
  const Vector< BaseIVFAB< Real>* >* m_coveredAdvVelPlusPtr;

  bool m_isNormalVelSet;
  bool m_isAdvVelSet;
  bool m_isMaxMinSet;
  Real m_maxVal;
  Real m_minVal;
private:

  BaseIVFAB<Vector<VolIndex> >  m_vofsStencil;

  //disallowed for all the usual reasons
  void operator=(const EBPatchTransport& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchTransport(const EBPatchTransport& a_input)
  {
    MayDay::Error("invalid operator");
  }

};


///
/**
 */
class EBPatchTransportFactory: public EBPatchGodunovFactory
{
public:
  ///
  /**
   */
  virtual ~EBPatchTransportFactory();

  ///
  /**
   */
  EBPatchGodunov* create() const;

  ///
  /**
   */
  EBPatchTransportFactory(RefCountedPtr<EBPhysIBCFactory>&  a_bc,
                       const bool&                       a_useLimiting,
                       const Real&                       a_maxVal =-1.e99,
                       const Real&                       a_minVal = 1.e99,
                       const bool&                       a_setMaxMin = false,
                       const bool&                       a_useFourthOrderSlopes = false,
                       const bool&                       a_useFlattening = false,
                       const bool&                       a_useArtificialVisc = false);

protected:
  RefCountedPtr<EBPhysIBCFactory>  m_bcFactoryPtr;
  bool                             m_useLimiting;
  bool                             m_setMaxMin;
  Real                             m_maxVal;
  Real                             m_minVal;
  bool                             m_useFourthOrderSlopes;
  bool                             m_useFlattening;
  bool                             m_useArtificialVisc;
  BaseIVFAB<Vector<VolIndex> >  m_vofsStencil;

private:


  //disallowed for all the usual reasons
  EBPatchTransportFactory()
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const EBPatchTransportFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchTransportFactory(const EBPatchTransportFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
