#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _KEPSILONPARAMS_H_
#define _KEPSILONPARAMS_H_

#include "NamespaceHeader.H"
#include "EBAMRRANSParams.H"

class kEpsilonParams : public EBAMRRANSParams
{
public:

  ~kEpsilonParams()
  {;}

  bool m_limitK;
  bool m_limitEpsilon;

  Real m_maxK;
  Real m_minK;
  Real m_minEpsilon;
  Real m_maxEpsilon;

  kEpsilonParams() : EBAMRRANSParams()
  {
    m_limitK = false;
    m_limitEpsilon = false;
    m_maxK = -1.e99;
    m_minK = 1.e99;
    m_maxEpsilon = -1.e99;
    m_minEpsilon = 1.e99;
  }
};

#include "NamespaceFooter.H"
#endif
