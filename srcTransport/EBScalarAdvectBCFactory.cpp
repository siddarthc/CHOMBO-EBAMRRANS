#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBScalarAdvectBCFactory.H"
#include "EBScalarAdvectBC.H"
#include "EBScalarBCF_F.H"

#include "NamespaceHeader.H"

/******************/
EBScalarAdvectBCFactory::
EBScalarAdvectBCFactory(int                                    a_inflowDir,
                        Real                                   a_scalarInflowValue,
                        RealVect                               a_injectOrigin,
                        RealVect                               a_injectEnd,
                        bool                                   a_initScalar,
                        RefCountedPtr<PoiseuilleInflowBCValue> a_injectionBCFunc)
  :EBPhysIBCFactory()
{
  m_inflowDir         = a_inflowDir;
  m_scalarInflowValue = a_scalarInflowValue;
  m_injectOrigin      = a_injectOrigin;  // for initial condition
  m_injectEnd         = a_injectEnd; // for initial condition
  m_injectionBCFunc   = a_injectionBCFunc;
  m_initScalar        = a_initScalar;

  Real injectRad = m_injectionBCFunc->getTubeRadius();
  Vector<Real> injectOrig(CH_SPACEDIM), injectEnd(CH_SPACEDIM);

  for (int i = 0; i < CH_SPACEDIM; i++)
    {
      injectOrig[i] = m_injectOrigin[i];
      injectEnd[i] = m_injectEnd[i];
    }

  if (m_initScalar)
    {
      FORT_SETINJECTIONPARAMS(CHF_CONST_INT(m_inflowDir),
                              CHF_CONST_REAL(m_scalarInflowValue),
                              CHF_CONST_REAL(injectRad),
                              CHF_CONST_VR(injectOrig),
                              CHF_CONST_VR(injectEnd));
      
      m_isFortranCommonSet = true;
    } 
}
/******************/
EBScalarAdvectBCFactory::
 ~EBScalarAdvectBCFactory()
{;}
/******************/
EBPhysIBC*
EBScalarAdvectBCFactory::
create() const
{
  if (m_initScalar) CH_assert(m_isFortranCommonSet);

  EBScalarAdvectBC* retval = new EBScalarAdvectBC(m_inflowDir, m_scalarInflowValue, m_injectOrigin, m_injectEnd, m_initScalar, m_injectionBCFunc);

  return static_cast<EBPhysIBC*>(retval);
}
/******************/

#include "NamespaceFooter.H"
