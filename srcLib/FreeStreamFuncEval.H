#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _FREESTREAMFUNCEVAL_H_
#define _FREESTREAMFUNCEVAL_H_

#include "BaseBCFuncEval.H"

#include "NamespaceHeader.H"

class FreeStreamFuncEval: public BaseBCFuncEval
{
public:
  
  FreeStreamFuncEval(const Real& a_freeStreamVel, const int& a_flowDir)
  {
    m_freeStreamVel = a_freeStreamVel;
    m_flowDir = a_flowDir;
  }

  ~FreeStreamFuncEval()
  {;}

  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp,
                     const Real&           a_time) const
  {
    Real retval = 0.;
    if (a_comp == m_flowDir) retval = m_freeStreamVel;
    return retval;
  }

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir) const
  {
    // hardwiring this
    return 0.;
  }

  Real m_freeStreamVel;
  int  m_flowDir;
};

#include "NamespaceFooter.H"
#endif

