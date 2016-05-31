#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBLevelRANS.H"
#include "EBPatchTransport.H"

#include "NamespaceHeader.H"

/***********/
EBLevelRANS::
EBLevelRANS() : EBLevelTransport()
{

}
/**********/
EBLevelRANS::
~EBLevelRANS()
{
}
/**********/
EBLevelRANS::
EBLevelRANS(const DisjointBoxLayout&           a_thisDBL,
            const DisjointBoxLayout&           a_coarDBL,
            const EBISLayout&                  a_thisEBISL,
            const EBISLayout&                  a_coarEBISL,
            const ProblemDomain&               a_domain,
            const int&                         a_nRefine,
            const RealVect&                    a_dx,
            const bool&                        a_hasCoarser,
            const bool&                        a_hasFiner,
            RefCountedPtr<EBPatchRANSModel>&   a_patchRANS,
            const EBIndexSpace*          const a_eb)
{
  CH_TIME("EBLevelRANS::EBLevelRANS");
  m_isDefined = false;
  define(a_thisDBL,
         a_coarDBL,
         a_thisEBISL,
         a_coarEBISL,
         a_domain,
         a_nRefine,
         a_dx,
         a_hasCoarser,
         a_hasFiner,
         a_patchRANS,
         a_eb);
}
/**********/
void 
EBLevelRANS::
define(const DisjointBoxLayout&           a_thisDBL,
       const DisjointBoxLayout&           a_coarDBL,
       const EBISLayout&                  a_thisEBISL,
       const EBISLayout&                  a_coarEBISL,
       const ProblemDomain&               a_DProblem,
       const int&                         a_nRefine,
       const RealVect&                    a_dx,
       const bool&                        a_hasCoarser,
       const bool&                        a_hasFiner,
       RefCountedPtr<EBPatchRANSModel>&   a_patchRANS,
       const EBIndexSpace*                const a_eb)
{
  RefCountedPtr<EBPatchTransport> patchTransport = static_cast<RefCountedPtr<EBPatchTransport> >(a_patchRANS);

  EBLevelTransport::define(a_thisDBL,a_coarDBL,a_thisEBISL,a_coarEBISL,a_DProblem,a_nRefine,a_dx,a_hasCoarser,a_hasFiner,patchTransport,a_eb);
}
/**********/
#include "NamespaceFooter.H"
