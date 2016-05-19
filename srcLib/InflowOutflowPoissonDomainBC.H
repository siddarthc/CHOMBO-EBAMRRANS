#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INFLOWOUTFLOWPOISSONDOMAINBC_H_
#define _INFLOWOUTFLOWPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "PoiseuilleInflowBCValue.H"

#include "NamespaceHeader.H"

///
/**
 */
class InflowOutflowPoissonDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~InflowOutflowPoissonDomainBC()
  {;}

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Implementation of the query of whether the face between a_ivof and a_jvof Dirichlet Domain boundary
     For InflowOutflowPoissonDomainBC, a_phi and a_dit are not used.
     For the moment it is only valid for the pressure calculation
     a_ivof must be inside the valid domain
  */
  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi) const;
  ///
  /**
     Elliptic solver inhom flux.
  */
  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
                                const VolIndex&       a_vof,
                                const int&            a_comp,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const int&            a_idir,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a InflowOutflow  value.
  */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);

  ///
  InflowOutflowPoissonDomainBC(const int& a_flowDir,
                               const bool& a_doJet2,
                               const Real& a_jet1inflowVel,
                               const Real& a_jet2inflowVel,
                               const bool& a_doJet1PoiseInflow,
                               const bool& a_doJet2PoiseInflow,
                               const IntVect& a_doSlipWallsHi,
                               const IntVect& a_doSlipWallsLo,
                               const ProblemDomain& a_domain,
                               const RealVect& a_dx,
                               RefCountedPtr<PoiseuilleInflowBCValue> a_jet1PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                               RefCountedPtr<PoiseuilleInflowBCValue> a_jet2PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> ())
  {
    m_flowDir = a_flowDir;
    m_doJet2  = a_doJet2;
    m_jet1inflowVel = a_jet1inflowVel;
    m_jet2inflowVel = a_jet2inflowVel;
    m_doJet1PoiseInflow = a_doJet1PoiseInflow;
    m_doJet2PoiseInflow = a_doJet2PoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_domain = a_domain;
    m_dx = a_dx;
    m_jet1PoiseInflowFunc = a_jet1PoiseBCValue;
    m_jet2PoiseInflowFunc = a_jet2PoiseBCValue;
  }

private:

  ///weak construction forbidden to keep things simple
  InflowOutflowPoissonDomainBC()
  {
  }

  int m_flowDir;
  bool m_doJet2;
  Real m_jet1inflowVel;
  Real m_jet2inflowVel;
  bool m_doJet1PoiseInflow;
  bool m_doJet2PoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  ProblemDomain m_domain;
  RealVect m_dx;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet1PoiseInflowFunc;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet2PoiseInflowFunc;

};

class InflowOutflowPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:


  ///
  /**
   */
  InflowOutflowPoissonDomainBCFactory(const int& a_flowDir,
                                      const bool& a_doJet2,
                                      const Real& a_jet1inflowVel,
                                      const Real& a_jet2inflowVel,
                                      const bool& a_doJet1PoiseInflow,
                                      const bool& a_doJet2PoiseInflow,
                                      const IntVect& a_doSlipWallsHi,
                                      const IntVect& a_doSlipWallsLo,
                                      RefCountedPtr<PoiseuilleInflowBCValue> a_jet1PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                                      RefCountedPtr<PoiseuilleInflowBCValue> a_jet2PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> ())
  {
    m_flowDir   = a_flowDir;
    m_doJet2    = a_doJet2;
    m_jet1inflowVel = a_jet1inflowVel;
    m_jet2inflowVel = a_jet2inflowVel;
    m_doJet1PoiseInflow = a_doJet1PoiseInflow;
    m_doJet2PoiseInflow = a_doJet2PoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_jet1PoiseInflowFunc = a_jet1PoiseBCValue;
    m_jet2PoiseInflowFunc = a_jet2PoiseBCValue;
  }


  ///
  /**
   */
  virtual ~InflowOutflowPoissonDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual InflowOutflowPoissonDomainBC* create(const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx)
  {
    InflowOutflowPoissonDomainBC* newBC
      = new InflowOutflowPoissonDomainBC(m_flowDir,
                                         m_doJet2,
                                         m_jet1inflowVel,
                                         m_jet2inflowVel,
                                         m_doJet1PoiseInflow,
                                         m_doJet2PoiseInflow,
                                         m_doSlipWallsHi,
                                         m_doSlipWallsLo,
                                         a_domain,
                                         a_dx,
                                         m_jet1PoiseInflowFunc,
                                         m_jet2PoiseInflowFunc);

    return newBC;
  }


private:
  InflowOutflowPoissonDomainBCFactory()
  {
  }

  int m_flowDir;
  bool m_doJet2;
  Real m_jet1inflowVel;
  Real m_jet2inflowVel;
  bool m_doJet1PoiseInflow;
  bool m_doJet2PoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet1PoiseInflowFunc;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet2PoiseInflowFunc;

};



///
/**
 */
class InflowOutflowHelmholtzDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~InflowOutflowHelmholtzDomainBC()
  {;}

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver inhom flux.
  */
  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
                                const VolIndex&       a_vof,
                                const int&            a_comp,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const int&            a_idir,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi) const;
  ///
  /**
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a InflowOutflow  value.
  */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);


  ///
  InflowOutflowHelmholtzDomainBC(const int& a_flowDir,
                                 const bool& a_doJet2,
                                 const Real& a_jet1inflowVel,
                                 const Real& a_jet2inflowVel,
                                 const bool& a_doJet1PoiseInflow,
                                 const bool& a_doJet2PoiseInflow,
                                 const IntVect& a_doSlipWallsHi,
                                 const IntVect& a_doSlipWallsLo,
                                 const ProblemDomain& a_domain,
                                 const RealVect& a_dx,
                                 RefCountedPtr<PoiseuilleInflowBCValue> a_jet1PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                                 RefCountedPtr<PoiseuilleInflowBCValue> a_jet2PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> ())
  {
    m_flowDir   = a_flowDir;
    m_doJet2    = a_doJet2;
    m_jet1inflowVel = a_jet1inflowVel;
    m_jet2inflowVel = a_jet2inflowVel;
    m_doJet1PoiseInflow = a_doJet1PoiseInflow;
    m_doJet2PoiseInflow = a_doJet2PoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_domain = a_domain;
    m_dx = a_dx;
    m_jet1PoiseInflowFunc = a_jet1PoiseBCValue;
    m_jet2PoiseInflowFunc = a_jet2PoiseBCValue;
  }

protected:
  //weak construction forbidden
  InflowOutflowHelmholtzDomainBC()
  {
    MayDay::Error("invalid operator");
  }

  int m_flowDir;
  bool m_doJet2;
  Real m_jet1inflowVel;
  Real m_jet2inflowVel;
  bool m_doJet1PoiseInflow;
  bool m_doJet2PoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  ProblemDomain m_domain;
  RealVect m_dx;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet1PoiseInflowFunc;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet2PoiseInflowFunc;
};

///
/**
 **/
class InflowOutflowHelmholtzDomainBCFactory: public BaseDomainBCFactory
{
public:

  ///
  InflowOutflowHelmholtzDomainBCFactory(const int& a_flowDir,
                                        const bool& a_doJet2,
                                        const Real& a_jet1inflowVel,
                                        const Real& a_jet2inflowVel,
                                        const bool& a_doJet1PoiseInflow,
                                        const bool& a_doJet2PoiseInflow,
                                        const IntVect& a_doSlipWallsHi,
                                        const IntVect& a_doSlipWallsLo,
                                        RefCountedPtr<PoiseuilleInflowBCValue> a_jet1PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                                        RefCountedPtr<PoiseuilleInflowBCValue> a_jet2PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> ())
  {
    m_flowDir   = a_flowDir;
    m_doJet2    = a_doJet2;
    m_jet1inflowVel = a_jet1inflowVel;
    m_jet2inflowVel = a_jet2inflowVel;
    m_doJet1PoiseInflow = a_doJet1PoiseInflow;
    m_doJet2PoiseInflow = a_doJet2PoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_jet1PoiseInflowFunc = a_jet1PoiseBCValue;
    m_jet2PoiseInflowFunc = a_jet2PoiseBCValue;
  }

  ///
  /**
   */
  virtual ~InflowOutflowHelmholtzDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual InflowOutflowHelmholtzDomainBC* create(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx)
  {
    InflowOutflowHelmholtzDomainBC* newBC
      = new InflowOutflowHelmholtzDomainBC(m_flowDir,
                                           m_doJet2,
                                           m_jet1inflowVel,
                                           m_jet2inflowVel,
                                           m_doJet1PoiseInflow,
                                           m_doJet2PoiseInflow,
                                           m_doSlipWallsHi,
                                           m_doSlipWallsLo,
                                           a_domain,
                                           a_dx,
                                           m_jet1PoiseInflowFunc,
                                           m_jet2PoiseInflowFunc);
    return newBC;
  }

protected:
  //weak construction
  InflowOutflowHelmholtzDomainBCFactory()
  {
    MayDay::Error("invalid operator");
  }

  int m_flowDir;
  bool m_doJet2;
  Real m_jet1inflowVel;
  Real m_jet2inflowVel;
  bool m_doJet1PoiseInflow;
  bool m_doJet2PoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet1PoiseInflowFunc;
  RefCountedPtr<PoiseuilleInflowBCValue> m_jet2PoiseInflowFunc;
};

#include "NamespaceFooter.H"

#endif
