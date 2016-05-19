#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBScalarAdvectBCFACTORY_H
#define _EBScalarAdvectBCFACTORY_H

#include "EBScalarAdvectBC.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

class EBScalarAdvectBCFactory: public EBPhysIBCFactory
{
public:

  ~EBScalarAdvectBCFactory();

  EBScalarAdvectBCFactory(int                                    a_inflowDir,
                          Real                                   a_scalarInflowValue,
                          RealVect                               a_injectOrigin,
                          RealVect                               a_injectEnd,
                          bool                                   a_initScalar,
                          RefCountedPtr<PoiseuilleInflowBCValue> a_injectionBCFunc = RefCountedPtr<PoiseuilleInflowBCValue> ());

  EBPhysIBC* create() const;

protected:
  int m_inflowDir;
  Real m_scalarInflowValue;
  RealVect m_injectOrigin;
  RealVect m_injectEnd;
  RefCountedPtr<PoiseuilleInflowBCValue> m_injectionBCFunc;
  bool m_initScalar;
  bool m_isFortranCommonSet;
private:
  EBScalarAdvectBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const EBScalarAdvectBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBScalarAdvectBCFactory(const EBScalarAdvectBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
