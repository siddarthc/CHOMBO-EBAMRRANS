#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PoiseuilleInflowBCValue.H"
#include "PolyGeom.H"

#include "NamespaceHeader.H"

/***/
RealVect
PoiseuilleInflowBCValue::
getVel(const Real& a_radius) const
{
  // CH_assert(SpaceDim==2);
  Real radsq = a_radius*a_radius;
  Real tubeRadSq = m_tubeRadius*m_tubeRadius;
  //    V = (1-r^2/r^2_0)*V_0
  Real vectSize = (1.0- (radsq/tubeRadSq))*m_maxVel;
  RealVect retval;
  // SC: hardwiring this 
  if (radsq <= tubeRadSq)
    {
      retval = m_tubeAxis*vectSize;
    }
  else 
    { 
      retval = RealVect::Zero;
    }
  return retval;
}
/***/
RealVect
PoiseuilleInflowBCValue::
getGradP() const
{
  CH_assert(SpaceDim==2);
  //    dp/dx = -8*[nu]/R^2 * u_avg
  Real tubeRadSq = m_tubeRadius*m_tubeRadius;
  Real vectSize = -(8.0*m_viscosity/(tubeRadSq))*m_avgVel;
  RealVect retval = m_tubeAxis*vectSize;
  return retval;
}
/***/
Real
PoiseuilleInflowBCValue::
getRadius(const RealVect& a_pt) const
{
  RealVect trueNorm, closestPt;
  PolyGeom::pointToLine(closestPt,  trueNorm, a_pt,
                        m_centerPt, m_tubeAxis);
  RealVect ebRadVec = closestPt;
  ebRadVec -= a_pt;

  Real ebRadius = PolyGeom::dot(ebRadVec, ebRadVec);
  ebRadius = sqrt(ebRadius);

  return ebRadius;
}
/***/
RealVect
PoiseuilleInflowBCValue::
getRadiusVector(const RealVect& a_pt) const
{
  RealVect trueNorm, closestPt;
  PolyGeom::pointToLine(closestPt,  trueNorm, a_pt,
                        m_centerPt, m_tubeAxis);
  RealVect ebRadVec = closestPt;
  ebRadVec -= a_pt;
  ebRadVec *= -1.;//wrong orientation

  return ebRadVec;
}
/***/
Real
PoiseuilleInflowBCValue::
value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const
{
  Real radius = this->getRadius(a_point);
  RealVect velocity = getVel(radius);
  return velocity[m_velComp];
}
/***/
Real
PoiseuilleInflowBCValue::
getViscosity() const
{
  return m_viscosity;
}
/***/
RealVect
PoiseuilleInflowBCValue::
getTubeCenter() const
{
  return m_centerPt;
}
/***/
Real 
PoiseuilleInflowBCValue::
getMaxVel() const
{
  return m_maxVel;
}
/***/
Real
PoiseuilleInflowBCValue::
getTubeRadius() const
{
  return m_tubeRadius;
}
#include "NamespaceFooter.H"
