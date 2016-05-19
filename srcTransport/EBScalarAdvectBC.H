#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EB_ScalarAdvectBC_H__
#define _EB_ScalarAdvectBC_H__

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "PoiseuilleInflowBCValue.H"
#include "EBPhysIBC.H"

#include "NamespaceHeader.H"

/**
  sets = inflow scalar value at inflow face.
  sets = zero at noflow faces.
  sets = extrapval at outflow face
*/

class EBScalarAdvectBC : public EBPhysIBC
{
public:
  virtual ~EBScalarAdvectBC();
  
  // specialized constructor for counter flow jet config
  EBScalarAdvectBC(int                                    a_inflowDir,
                   Real                                   a_scalarInflowValue,
                   RealVect                               a_injectOrigin,
                   RealVect                               a_injectEnd,
                   bool                                   a_initScalar,
                   RefCountedPtr<PoiseuilleInflowBCValue> a_injectionBCFunc = RefCountedPtr<PoiseuilleInflowBCValue> ());

  ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

  ///  For every box in this level, this function is called.
  void fluxBC(EBFluxFAB&            a_flux,
              const EBCellFAB&      a_Wcenter,
              const EBCellFAB&      a_Wextrap,
              const Side::LoHiSide& a_side,
              const Real&           a_time,
              const EBISBox&        a_ebisBox,
              const DataIndex&      a_dit,
              const Box&            a_box,
              const Box&            a_faceBox,
              const int&            a_faceDir);

  /// Initialize
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  void setAdvVel(const EBFluxFAB& a_advVel);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
  */
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir)
  {;}

protected:
  bool m_isDefined;
  Real m_inflowDir;
  Real m_scalarInflowValue;
  RealVect m_injectOrigin;
  RealVect m_injectEnd;
  bool m_initScalar;
  RefCountedPtr<PoiseuilleInflowBCValue> m_injectionBCFunc;

  RealVect m_dx;
  ProblemDomain  m_domain;

  bool m_isAdvVelSet;
  const EBFluxFAB* m_advVelPtr;
  
private:
  EBScalarAdvectBC()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBScalarAdvectBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBScalarAdvectBC(const EBScalarAdvectBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
