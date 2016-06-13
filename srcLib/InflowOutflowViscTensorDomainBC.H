#ifdef CH_LANG_CC
/*
 *        _______              __
 *       / ___/ /  ___  __ _  / /  ___
 *      / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *      \___/_//_/\___/_/_/_/_.__/\___/
 *      Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _InflowOutflowViscTensorDomainBC_H_
#define _InflowOutflowViscTensorDomainBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "PoiseuilleInflowBCValue.H"

#include "NamespaceHeader.H"

class InflowOutflowViscTensorDomainBC: public ViscousBaseDomainBC
{
public:

  virtual ~InflowOutflowViscTensorDomainBC()
  {;}

  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

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
                              const bool&           a_useHomogeneous)
  {
    //this is for projections and such.
    //not used here.
    MayDay::Error();
  }

  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow)
  {
    //this is for projections and such.
    //not used here.
    MayDay::Error();
  }

  InflowOutflowViscTensorDomainBC(const int& a_flowDir,
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
  InflowOutflowViscTensorDomainBC()
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
  
class InflowOutflowViscTensorDomainBCFactory: public BaseDomainBCFactory
{
public:

  ///
  InflowOutflowViscTensorDomainBCFactory(const int& a_flowDir,
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
  virtual ~InflowOutflowViscTensorDomainBCFactory()
  {
  }

  virtual InflowOutflowViscTensorDomainBC* create(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx)
  {
    InflowOutflowViscTensorDomainBC* newBC
      = new InflowOutflowViscTensorDomainBC(m_flowDir,
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
  InflowOutflowViscTensorDomainBCFactory()
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
