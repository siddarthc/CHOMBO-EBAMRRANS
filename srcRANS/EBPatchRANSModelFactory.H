#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHRANSMODELFACTORY_H_
#define _EBPATCHRANSMODELFACTORY_H_

#include "EBPatchRANSModel.H"
#include "EBAMRRANSParams.H"
#include "NamespaceHeader.H"

class EBPatchRANSModelFactory
{
public:

  virtual ~EBPatchRANSModelFactory();

  EBPatchRANSModel* create() const;

  EBPatchRANSModelFactory(RefCountedPtr<EBPhysIBCFactory>& a_bc,
                          EBAMRRANSParams& a_params,
                          const bool&                      a_useLimiting,
                          const bool&                      a_useFourthOrderSlopes,
                          const bool&                      a_useFlattening,
                          const bool&                      a_useArtificialVisc );

protected:

  RefCountedPtr<EBPhysIBCFactory>  m_bcFactoryPtr;
  EBAMRRANSParams m_params;
  bool m_useFourthOrderSlopes;
  bool m_useFlattening;
  bool m_useArtificialVisc;
  bool m_useLimiting;

private:

  EBPatchRANSModelFactory()
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const EBPatchRANSModelFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchRANSModelFactory(const EBPatchRANSModelFactory& a_input)
  {
    MayDay::Error("invalid operator");
  } 
};
#include "NamespaceFooter.H"
#endif
