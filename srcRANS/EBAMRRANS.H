#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRRANS_H_
#define _EBAMRRANS_H_

#include "EBAMRLevel.H"
#include "EBAMRRANSParams.H"
#include "EBPatchRANSModelFactory.H"
#include "EBPatchRANSModel.H"
#include "EBLevelRANS"

#include "NamespaceHeader.H"

// RANS solver 
// driven by external AMR or EBAMR class
class EBAMRRANS : public EBAMRLevel
{
public:

  virtual ~EBAMRRANS();

  EBAMRRANS(const EBAMRRANSParams& a_params,
            const RefCountedPtr<EBPatchRANSModelFactory>& a_modelFactory,
            bool  a_externalDriver);

};

#include "NamespaceFooter.H"
#endif
