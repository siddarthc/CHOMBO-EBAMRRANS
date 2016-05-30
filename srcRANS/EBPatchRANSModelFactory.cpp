#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchRANSModel.H"
#include "kEpsilon.H"
#include "kEpsilonParams.H"
#include "EBPatchRANSModelFactory.H"

#include "NamespaceHeader.H"

/*****/
EBPatchRANSModelFactory::
EBPatchRANSModelFactory(RefCountedPtr<EBPhysIBCFactory>& a_bcFactoryPtr,
                        EBAMRRANSParams&                 a_params,
                        const bool&                      a_useLimiting,
                        const bool&                      a_useFourthOrderSlopes,
                        const bool&                      a_useFlattening,
                        const bool&                      a_useArtificialVisc)
{
  m_bcFactoryPtr = a_bcFactoryPtr;
  m_params = a_params;
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_useFlattening = a_useFlattening;
  m_useArtificialVisc = a_useArtificialVisc;
  m_useLimiting = a_useLimiting;
}
/******/
EBPatchRANSModelFactory::
~EBPatchRANSModelFactory()
{
}
/*****/
EBPatchRANSModel*
EBPatchRANSModelFactory::
create() const
{
  // specializing for kEpsilon model
  // this is where the instantiantiation of appropriate model depending on 
  //                                the user input should take place 

  kEpsilon* retval = new kEpsilon();
  retval->setEBPhysIBC(*m_bcFactoryPtr);
  const kEpsilonParams* kEpsParams = dynamic_cast<const kEpsilonParams*>(&m_params);
  retval->setParams(*kEpsParams);
  retval->setSlopeParameters(m_useFourthOrderSlopes,
                             m_useFlattening,
                             m_useLimiting);
  retval->artificialViscosity(m_useArtificialVisc);
  return static_cast<EBPatchRANSModel*>(retval); 
}
/***************/
#include "NamespaceFooter.H"
