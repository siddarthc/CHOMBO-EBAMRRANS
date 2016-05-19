#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBAMRLevel.H"
#include "EBAMRTransport.H"
#include "EBAMRTransportFactory.H"

#include "NamespaceHeader.H"

/****************/
EBAMRLevel*
EBAMRTransportFactory::
new_amrlevel() const
{
  EBAMRTransport* amrg_ptr = new EBAMRTransport(m_params, m_patchFact, m_externalDriver);
  return (static_cast <EBAMRLevel*> (amrg_ptr));
}
/****************/

#include "NamespaceFooter.H"
