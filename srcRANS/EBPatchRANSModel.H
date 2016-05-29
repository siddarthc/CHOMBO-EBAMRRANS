#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHRANSMODEL_H_
#define _EBPATCHRANSMODEL_H_

#include "EBPatchTransport.H"
#include "NamespaceHeader.H"

#define TOL 1.e-12

///
/**
  Base class for the RANS model
  This essentially drives the individual RANS transport equations (k and epsilon for k-epsilon model, etc.)
*/
class EBPatchRANSModel : public EBPatchTransport
{
public:

  ///
  EBPatchRANSModel();

  ///
  ~EBPatchRANSModel();

  ///
  /**
    Return number of equations in the Model (2 for k-epsilpon)
  */
  virtual int numEquations() const = 0;  

  /**
     Return number of components for primitive variables.
  */
  virtual int numPrimitives() const
  {
    MayDay::Error("numPrimitives() should not be called from EBPatchRANSModel");
    return 0;
  }

  ///
  /**
     Returns number of components for flux variables.
  */
  virtual int numFluxes() const
  {
    MayDay::Error("numFluxes() should not be called from EBPatchRANSModel");
    return 0;
  }

  ///
  /**
     Returns number of components for flux variables.
  */
  virtual int numSlopes() const
  {
    MayDay::Error("numSlopes() should not be called from EBPatchRANSModel");
    return 0;
  }

  ///
  /**
     Returns number of components for conserved variables.
  */
  virtual int numConserved() const
  {
    MayDay::Error("numConserved() should not be called from EBPatchRANSModel");
    return 0;
  }

  ///
  /**
   */
  virtual void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box);   

   /// floors if m_isMaxMinSet 
   virtual void floorPrimitives(EBCellFAB& a_primState,
                                const Box& a_box)
   {
     MayDay::Error("floorPrimitives should not be called from EBPatchRANSModel");
   }

   /// floors if m_isMaxMinSet 
   virtual void floorPrimitives(BaseIVFAB<Real>&  a_primState,
                                const IntVectSet& a_ivsIrreg)
   {
     MayDay::Error("floorPrimitives should not be called from EBPatchRANSModel");
   }

  /// floors if m_isMaxMinSet
  virtual void
  floorConserved(EBCellFAB&  a_consState,
                 const Box&  a_box)
  {
     MayDay::Error("floorConserved should not be called from EBPatchRANSModel");
  }

  /// floors if m_isMaxMinSet
  virtual void
  floorConserved(BaseIVFAB<Real>&   a_consState,
                 const IntVectSet&  a_ivsIrreg)
  {
     MayDay::Error("floorConserved should not be called from EBPatchRANSModel");
  }

  virtual Vector<string> stateNames()
  {
     MayDay::Error("stateNames() should not be called from EBPatchRANSModel");
  }

  virtual Vector<string> primNames()
  {
     MayDay::Error("primNames() should not be called from EBPatchRANSModel");
  }

  virtual void
  getFlux(EBFluxFAB&       a_flux,
          const EBFluxFAB& a_prim);

  virtual void setMaxMin(const Real& a_maxVal,
                         const Real& a_minVal)
  {
     MayDay::Error("setMaxMin should not be called from EBPatchRANSModel");
  }

    // The essence 
  virtual void getTurbulentViscosity(EBFluxFAB&       a_nuT,
                                     const EBCellFAB& a_state) = 0;

protected:

  /// 
  /** 
    Set number of equations in the Model (2 for k-epsilpon, etc.)
  */
  virtual void setNumEquations() = 0;

  // number of equations
  int m_nEqn;

  int m_nPrim;
  int m_nCons;
  int m_nFlux;
  int m_nSlope;

  Real m_maxNuT;  //NuT being the turbulent viscosity
  Real m_minNuT;
  bool m_limitNuT; 

private:

  //disallowed for all the usual reasons
  void operator=(const EBPatchRANSModel& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchRANSModel(const EBPatchRANSModel& a_input)
  {
    MayDay::Error("invalid operator");
  }
}; 

#include "NamespaceFooter.H"
#endif 
