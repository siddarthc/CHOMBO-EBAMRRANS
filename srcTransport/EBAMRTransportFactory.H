#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRTRANSPORTFACTORY_H_
#define _EBAMRTRANSPORTFACTORY_H_

#include "EBAMRTransport.H"
#include "EBAMRTransportParams.H"
#include "EBAMRLevelFactory.H"
#include "EBPatchTransport.H"
#include "EBPatchGodunovFactory.H"

#include "NamespaceHeader.H"

class EBAMRTransportFactory: public EBAMRLevelFactory
{
  
  public:
  //! 
  EBAMRLevel* new_amrlevel() const;

  //! Create a new factory that assembles EBAMRTransport level solvers.
  //! \param a_params The parameters for the level solvers.
  //! \param a_factory A factory that builds EBPatchTransport objects for the solvers.
  EBAMRTransportFactory(const EBAMRTransportParams& a_params,
                        const RefCountedPtr<EBPatchGodunovFactory>& a_patchFact,
                        bool a_extranalDriver):
    m_params(a_params),
    m_patchFact(a_patchFact),
    m_externalDriver(a_extranalDriver)
  {
  }

  virtual ~EBAMRTransportFactory()
  {;}

  protected:

  EBAMRTransportParams m_params;
  RefCountedPtr<EBPatchGodunovFactory> m_patchFact;
  bool m_externalDriver;

  private:
  
  // Disallowed operators.
  EBAMRTransportFactory();

  // The compiler generates the copy constructor and assignment operator. 
};

#include "NamespaceFooter.H"
#endif
