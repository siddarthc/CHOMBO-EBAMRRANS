#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRTransportParams_H_
#define _EBAMRTransportParams_H_

#include "NamespaceHeader.H"

class EBAMRTransportParams
{
public:
  bool             m_doDiffusion;
  int              m_verbosity;
  int              m_redistRad;
  bool             m_useMassRedist;
  bool             m_doSmushing;
  bool             m_variableCoeff;
  bool             m_checkMaxMin;
  Real             m_diffusionCoeff; // this is the beta for conductivityOp
  Real             m_alphaCoeff; // this the the alpha for conductivityOp (rho*Cv in case of temperature)
  bool             m_doViscDissipation; 
  bool             m_backwardEuler;

  EBAMRTransportParams()
  {
    m_doDiffusion   = true;
    m_verbosity     = 3;
    m_redistRad     = 1;
    m_useMassRedist = true; // this is to set weights on the FinetoCoar and CoarToCoar objects based on the scalar value
    m_doSmushing    = false; // 
    m_variableCoeff = false;
    m_checkMaxMin   = true;
    m_alphaCoeff    = 1.0;
    m_backwardEuler = false;
    m_diffusionCoeff = 1.0;
    m_doViscDissipation = false;
  }

  ~EBAMRTransportParams()
  {;}
};

#include "NamespaceFooter.H"
#endif
