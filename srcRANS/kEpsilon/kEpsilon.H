#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _KEPSILON_H_
#define _KEPSILON_H_

#include "EBPatchRANSModel.H"
#include "kEpsilonParams.H"
#include "NamespaceHeader.H"

#define TOL 1.e-12

///
/**
  Base kEpsilon model
  This lives only on patches
  The other versions of kEpslion model may be derived from this
*/
class kEpsilon : public EBPatchRANSModel
{
public:

  ///
  kEpsilon();

  void setParams(const kEpsilonParams& a_params)
  {
    m_params = a_params;
  }

  ///
  ~kEpsilon();

  virtual int numEquations() const
  {
    return m_nEqn;
  }

  virtual int numPrimitives() const 
  {
    return m_nPrim; 
  } 

  virtual int numConserved() const 
  {
    return m_nCons;
  }

  virtual int numFluxes() const
  {
    return m_nFlux;
  }

  virtual int numSlopes() const
  {
    return m_nSlope;
  }

  /// floors if m_isMaxMinSet 
  virtual void floorPrimitives(EBCellFAB& a_primState,
                               const Box& a_box); 

  virtual void floorPrimitives(BaseIVFAB<Real>&  a_primState,
                                const IntVectSet& a_ivsIrreg);

  virtual void floorConserved(EBCellFAB&  a_consState,
                              const Box&  a_box);

  virtual void floorConserved(BaseIVFAB<Real>&   a_consState,
                              const IntVectSet&  a_ivsIrreg);
 
  virtual Vector<string> stateNames();

  virtual Vector<string> primNames(); 

  virtual void getTurbulentViscosity(EBFluxFAB&       a_nuT,
                                     const EBCellFAB& a_state)
  {
    // follow EBPatchTransport::primitivesAndDivergences
  }

  virtual void computeNetSource(EBCellFAB&       a_netSource,
                                const EBCellFAB& a_state,
                                const Box&       a_box);

  /// the equation index depending on which stencil weights are reset (look EBAMRRANS::resetWeights)
  virtual Real stencilWeightEqnIndex();

  virtual void getDiffusionCoefficients(EBFluxFAB&       a_diffCoeff,
                                        const EBFluxFAB& a_stateFace,
                                        const Box&       a_box);

  virtual void getDiffusionCoefficients(BaseIVFAB<Real>& a_diffCoeffIrreg,
                                        const EBCellFAB& a_stateCell,
                                        const IntVectSet& a_ivs);

    virtual void fillDriverDiffusionCoefficients(EBFluxFAB&       a_diffCoeff,
                                        const EBFluxFAB& a_stateFace,
                                        const Box&       a_box);

  virtual void fillDriverDiffusionCoefficients(BaseIVFAB<Real>& a_diffCoeffIrreg,
                                        const EBCellFAB& a_stateCell,
                                        const IntVectSet& a_ivs);

protected:
  
  /// 
  /**
    Set number of equations
  */
  virtual void setNumEquations();

  
  int m_kEqnIndex;  
  int m_epsEqnIndex;  
/*
  Real m_maxK;
  Real m_maxEpsilon;

  Real m_minK;
  Real m_minEpsilon;

  bool m_limitK;
  bool m_limitEpsilon;
*/
  kEpsilonParams m_params;
private:

  //disallowed for all the usual reasons
  void operator=(const kEpsilon& a_input)
  {
    MayDay::Error("invalid operator");
  }
  kEpsilon(const kEpsilon& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
