#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "Box.H"
#include "Vector.H"
#include "LayoutIterator.H"
#include "parstream.H"

#include "EBAMRLevel.H"
#include "NamespaceHeader.H"

bool EBAMRLevel::s_noEBCF = false;
int EBAMRLevel::s_verbosity = 0;

//-----------------------------------------------------------------------
bool EBAMRLevel::isDefined() const
{
  return m_isDefined;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBAMRLevel::~EBAMRLevel()
{
//  if (m_dtPtr != NULL) delete m_dtPtr;
//  if (m_timePtr != NULL) delete m_timePtr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBAMRLevel::EBAMRLevel()
{
  m_coarser_level_ptr = NULL;
  m_finer_level_ptr = NULL;
  m_isDefined = false;
  m_level = 0;
  m_time = 0;
  m_dt = 0;
  m_initial_dt_multiplier = 0.1;
//  m_isEBLGSet = false;
  m_oldEBLGSet = false;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::define(EBAMRLevel*  a_coarser_level_ptr,
                        const Box&   a_problem_domain,
                        int          a_level,
                        int          a_ref_ratio)
{
  ProblemDomain physDomain(a_problem_domain);

  define(a_coarser_level_ptr, physDomain, a_level, a_ref_ratio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::define(EBAMRLevel*            a_coarser_level_ptr,
                        const ProblemDomain&   a_problem_domain,
                        int                    a_level,
                        int                    a_ref_ratio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::define" << endl;
  }

  m_coarser_level_ptr = a_coarser_level_ptr;
  m_problem_domain = a_problem_domain;
  m_level = a_level;
  m_ref_ratio = a_ref_ratio;
  m_finer_level_ptr = NULL;
  m_isDefined = true;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::finerLevelPtr(EBAMRLevel* a_finer_level_ptr)
{
  m_finer_level_ptr = a_finer_level_ptr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::dt(Real a_dt)
{
  m_dt = a_dt;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::dtPtr(const Real* a_dtPtr)
{
  m_dtPtr = a_dtPtr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::timePtr(const Real* a_timePtr)
{
  m_timePtr = a_timePtr;
} 
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real EBAMRLevel::dt() const
{
  return(m_dt);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
const ProblemDomain& EBAMRLevel::problemDomain() const
{
  return(m_problem_domain);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<Box> EBAMRLevel::boxes() const
{
  return(m_level_grids);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool EBAMRLevel::hasCoarserLevel() const
{
  return ((m_coarser_level_ptr != NULL)
       && (m_coarser_level_ptr->m_isDefined)
       && (m_coarser_level_ptr->m_level_grids.size() > 0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool EBAMRLevel::hasFinerLevel() const
{
  return ((m_finer_level_ptr != NULL)
       && (m_finer_level_ptr->m_isDefined)
       && (m_finer_level_ptr->m_level_grids.size() > 0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int EBAMRLevel::level() const
{
  return m_level;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int EBAMRLevel::refRatio() const
{
  return(m_ref_ratio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::time(Real a_time)
{
  m_time = a_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real EBAMRLevel::time() const
{
  return m_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::initialDtMultiplier(Real a_initial_dt_multiplier)
{
  m_initial_dt_multiplier = a_initial_dt_multiplier;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real EBAMRLevel::initialDtMultiplier() const
{
  return m_initial_dt_multiplier;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// static
void EBAMRLevel::verbosity(int a_verbosity)
{
  s_verbosity = a_verbosity;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// static
int EBAMRLevel::verbosity()
{
  return(s_verbosity);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::preRegrid" << endl;
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::postRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::postRegrid" << endl;
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::postInitialGrid(const bool a_restart)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::postInitialGrid" << endl;
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::setEBLGOld(const EBLevelGrid& a_eblg)
{
  m_eblgOld = a_eblg;
  m_oldEBLGSet = true;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void EBAMRLevel::setTimeOld(Real& a_dt, Real& a_time)
{
  m_dtOld = a_dt;
  m_timeOld = a_time; 
  m_oldTimeSet = true;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
/*
void EBAMRLevel::setEBLG(const EBLevelGrid& a_eblg)
{
  m_eblg = a_eblg;
  m_isEBLGSet = true;
}
*/
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<EBAMRLevel*> EBAMRLevel::getAMRLevelHierarchy()
{
  Vector<EBAMRLevel*> retval;
  // First go to level 0
  EBAMRLevel* levelPtr = this;
  while (levelPtr->hasCoarserLevel())
    {
      levelPtr = levelPtr->m_coarser_level_ptr;
    }

  // Now can accumulate the pointers by chasing finer level
  retval.push_back(levelPtr);
  while (levelPtr->hasFinerLevel())
    {
      levelPtr = levelPtr->m_finer_level_ptr;
      retval.push_back(levelPtr);
    }

  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBAMRLevel::writeCustomPlotFile(const std::string& a_prefix,
                                int a_step) const
{
  // By default, this does nothing.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBAMRLevel::conclude(int a_step) const
{
  // By default, this does nothing.
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
