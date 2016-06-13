#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _POISEUILLEINFLOWFUNCEVAL_H_
#define _POISEUILLEINFLOWFUNCEVAL_H_

#include "BaseBCFuncEval.H"
#include "PoiseuilleInflowBCValue.H"

#include "NamespaceHeader.H"

// just a wrapper for viscous tensor op
class PoiseuilleInflowFuncEval: public BaseBCFuncEval
{
public:

  PoiseuilleInflowFuncEval(RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> ())
  {
    m_func = a_poiseBCValue;
  }

  ~PoiseuilleInflowFuncEval()
  {
  }

  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp,
                     const Real&           a_time) const
{
  Real radius = m_func->getRadius(a_point);
  RealVect velocity = m_func->getVel(radius);
  return velocity[a_comp]; 
}

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir) const
{
  // hardwiring this
  return 0.;
}

  RefCountedPtr<PoiseuilleInflowBCValue> m_func;

};

#include "NamespaceFooter.H"
#endif

