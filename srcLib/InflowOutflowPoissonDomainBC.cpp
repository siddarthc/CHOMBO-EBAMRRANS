#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "InflowOutflowPoissonDomainBC.H"
#include "PoiseuilleInflowBCValue.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletPoissonEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

extern Real g_simulationTime;
bool s_higherOrderHelmBC = false;

//pressure stencil
void InflowOutflowPoissonDomainBC::getFluxStencil(      VoFStencil&      a_stencil,
                                                  const VolIndex&        a_vof,
                                                  const int&             a_comp,
                                                  const RealVect&        a_dx,
                                                  const int&             a_idir,
                                                  const Side::LoHiSide&  a_side,
                                                  const EBISBox&         a_ebisBox)
{
  //added: 04/20 04:51
  bool isInsideTube = false;
  if ((m_doJet2) && (a_side==Side::Hi) && (a_idir == m_flowDir))
    {
      RealVect probLo = RealVect::Zero;
      RealVect loc = EBArith::getVofLocation(a_vof,a_dx,probLo);
      loc[a_idir] += 0.5*a_dx[a_idir]; //loc at face center 
      const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
      Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
      Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
      CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);

// start hardwire
      Real wallThickness = 0.075;
      Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
      ParmParse pp;
      pp.get("wall_thickness",wallThickness); 
      tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

      if (radius <= tubeRadius)
        {
          isInsideTube = true;
        }
    }
  bool isOutflow = (a_side==Side::Hi) && (a_idir==m_flowDir) && (!(isInsideTube));
  // end addition
//  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setEBOrder(1);
      diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
}

///
/**
   This is called by cc projection to enforceVelocityBCs on face
   velocities averaged from centers before calculating divergence in mac projection.
*/
void
InflowOutflowPoissonDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  CH_assert(a_vel.nComp() == 1);
  bool isInflow = (a_side==Side::Lo) && (a_idir==m_flowDir);
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);

  RealVect point =  EBArith::getFaceLocation(a_face, a_dx, a_probLo);
  RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

  const EBISBox& ebisBox= a_vel[0].getEBISBox();

  if (isOutflow)
    {
      if (!(m_doJet2))
        {
          //quadratic extrapolation with homogeneous Neumann for outflow bc
          a_faceFlux = EBArith::extrapFaceVelToOutflow(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0);//be careful of this 0, it is the component
        }
      else if (m_doJet2)
        {
          const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
          Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
          CH_assert(point[m_flowDir] == tubeCenter[m_flowDir]);
          Real radius = m_jet2PoiseInflowFunc->getRadius(point);
          bool isInsideTube = false;

// start hardwire
          Real wallThickness = 0.075;
          Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
          ParmParse pp;
          pp.get("wall_thickness",wallThickness);
          tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

          if (radius <= tubeRadius)
            {
              isInsideTube = true;
            }
          if (isInsideTube && (velcomp == m_flowDir))
            {
              if (m_doJet2PoiseInflow)
                {
                  a_faceFlux = m_jet2PoiseInflowFunc->getVel(radius)[m_flowDir];
                }
              else if (!(m_doJet2PoiseInflow))
                {
                  a_faceFlux = m_jet2inflowVel;
                }
            }
          else if (isInsideTube && (velcomp != m_flowDir))
            {
              a_faceFlux = 0; // set tangential components to 0
            }
          else if (!(isInsideTube))
            {
              //quadratic extrapolation with homogeneous Neumann for outflow bc
              a_faceFlux = EBArith::extrapFaceVelToOutflow(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0);//be careful of this 0, it is the component 
            } 
        }
    }
  else if (isInflow && (velcomp==m_flowDir))
    {
      //input component is always zero.
      //value of face direction corresponds to velocity direction here
      if (!m_doJet1PoiseInflow)
      {
        a_faceFlux = m_jet1inflowVel;
      }
      else if (m_doJet1PoiseInflow)
        {
          RealVect prob_lo = RealVect::Zero;
          const RealVect loc  = EBArith::getFaceLocation(a_face, a_dx, prob_lo);
          Real radius = m_jet1PoiseInflowFunc->getRadius(loc);
          a_faceFlux = m_jet1PoiseInflowFunc->getVel(radius)[m_flowDir];
        }
    }
  else if (isInflow && (velcomp != m_flowDir))
    {
      //tangential velocity at inflow
      a_faceFlux = 0;
    }
  else
    {
      // solid wall = faceFlux = 0; free stream = extrap (0 Neum)
      ParmParse pp;
      int bcType; // 0 solid wall, 1 free stream
      pp.get("bcType", bcType);
      if (bcType == 0)  // solid wall
       {
         a_faceFlux = 0;
       }
     else if (bcType == 1) // free stream
       {
         //quadratic extrapolation with homogeneous Neumann 
         a_faceFlux = EBArith::extrapFaceVelToOutflow(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0);//be careful of this 0, it is the component 
       }
     else
       {
         MayDay::Error("wrong BC type");
       } 
    }
}

//////called by MAC projection solver for domain bc on phi
void InflowOutflowPoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                               const BaseFab<Real>&  a_phi,
                                               const RealVect&       a_probLo,
                                               const RealVect&       a_dx,
                                               const int&            a_idir,
                                               const Side::LoHiSide& a_side,
                                               const DataIndex&      a_dit,
                                               const Real&           a_time,
                                               const bool&           a_useHomogeneous)
{

  //for phi: outflow  dirichlet. inflow neumann
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
  else
    {
      if (!(m_doJet2))
        {
          DirichletPoissonDomainBC diriBC;
          diriBC.setValue(0.0);
          diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
          // diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else if (m_doJet2)
        {
          const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
          Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();

// start hardwire
          Real wallThickness = 0.075;
          Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
          ParmParse pp;
          pp.get("wall_thickness",wallThickness);
          tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

          CH_assert(a_phi.nComp() == 1);
          for (int comp=0; comp<a_phi.nComp(); comp++)
            {
              const Box& box = a_faceFlux.box();

              int iside = -1; // a_side == Side::Hi in this loop
              Real value = 0.0;
              BoxIterator bit(box);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  RealVect loc = EBArith::getIVLocation(iv,a_dx,a_probLo);
                  loc[a_idir] -= iside * 0.5 * a_dx[a_idir]; //loc is now at the face center
                  CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);
                  bool isInsideTube = false;
                  Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
                  if (radius <= tubeRadius)
                    {
                      isInsideTube = true;
                    }
                   if (isInsideTube)
                    {
                      a_faceFlux(iv,comp) = iside * value; //Neumann
                    }
                  else if (!(isInsideTube))
                    {
                      Real ihdx = 2.0 / a_dx[a_idir];
                      Real phiVal = a_phi(iv,comp);
                      a_faceFlux(iv,comp) = iside * ihdx * (phiVal - value); //Dirichlet
                    }
                }
            }
        }
    }
}

//////called by MAC projection solver for domain bc on phi
bool InflowOutflowPoissonDomainBC::
isDirichletDom(const VolIndex&   a_ivof,
               const VolIndex&   a_jvof,
               const EBCellFAB&  a_phi) const
{

/*
  // we don't really need the FaceIndex object, just trying to save some codes here
  FaceIndex face(a_ivof,a_jvof);
  int side = face.faceSign(a_ivof);
  int idir = face.direction();

  bool isOutflow = (side==1) && (idir==m_flowDir);
  // this logic is valid only for the projection
  if (!isOutflow)
    {
      return false;
    }
  else
    {
      return true;
    }
*/
  MayDay::Error("InfowOutflowHelmholtzDomain::isDirichletDom is not set");
  return true; // to shut the compiler
}

//called by projection solve for irreg x domain
void InflowOutflowPoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
                                               const VolIndex&       a_vof,
                                               const int&            a_comp,
                                               const EBCellFAB&      a_phi,
                                               const RealVect&       a_probLo,
                                               const RealVect&       a_dx,
                                               const int&            a_idir,
                                               const Side::LoHiSide& a_side,
                                               const DataIndex&      a_dit,
                                               const Real&           a_time,
                                               const bool&           a_useHomogeneous)
{
  bool isInsideTube = false;
  if ((m_doJet2) && (a_side==Side::Hi) && (a_idir == m_flowDir))
    {
      RealVect probLo = RealVect::Zero;
      RealVect loc = EBArith::getVofLocation(a_vof,a_dx,probLo);
      loc[a_idir] += 0.5*a_dx[a_idir]; // loc at face center
      const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
      Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
      Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
      CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);

// start hardwire
      Real wallThickness = 0.075;
      Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
      ParmParse pp;
      pp.get("wall_thickness",wallThickness);
      tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

      if (radius <= tubeRadius)
        {
          isInsideTube = true;
        }
    }
  bool isOutflow = (a_side==Side::Hi) && (a_idir==m_flowDir) && (!(isInsideTube));

//  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      // diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
      //                               a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
}

//called by projection solve for irreg x domain
void InflowOutflowPoissonDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
                                                    const VolIndex&       a_vof,
                                                    const int&            a_comp,
                                                    const EBCellFAB&      a_phi,
                                                    const RealVect&       a_probLo,
                                                    const RealVect&       a_dx,
                                                    const int&            a_idir,
                                                    const Side::LoHiSide& a_side,
                                                    const DataIndex&      a_dit,
                                                    const Real&           a_time)
{
  bool isInsideTube = false;
  if ((m_doJet2) && (a_side==Side::Hi) && (a_idir == m_flowDir))
    {
      RealVect probLo = RealVect::Zero;
      RealVect loc = EBArith::getVofLocation(a_vof,a_dx,probLo);
      loc[a_idir] += 0.5*a_dx[a_idir]; // loc at face center
      const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
      Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
      Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
      CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);

// start hardwire
      Real wallThickness = 0.075;
      Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
      ParmParse pp;
      pp.get("wall_thickness",wallThickness);
      tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

      if (radius <= tubeRadius)
        {
          isInsideTube = true;
        }
    }
  bool isOutflow = (a_side==Side::Hi) && (a_idir==m_flowDir) && (!(isInsideTube));
//  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit, a_time);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                              a_dx, a_idir, a_side, a_dit, a_time);
      // diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
      //                               a_dx, a_idir, a_side, a_dit, a_time);
    }
}

//called by macEnforceGradientBC, useAreaFrac = false
void InflowOutflowPoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                  const FaceIndex&      a_face,
                                                  const int&            a_comp,
                                                  const EBCellFAB&      a_phi,
                                                  const RealVect&       a_probLo,
                                                  const RealVect&       a_dx,
                                                  const int&            a_idir,
                                                  const Side::LoHiSide& a_side,
                                                  const DataIndex&      a_dit,
                                                  const Real&           a_time,
                                                  const bool&           a_useAreaFrac,
                                                  const RealVect&       a_centroid,
                                                  const bool&           a_useHomogeneous)
{
  bool isInsideTube = false;
  if ((m_doJet2) && (a_side==Side::Hi) && (a_idir == m_flowDir))
    {
      RealVect probLo = RealVect::Zero;
      RealVect loc = EBArith::getFaceLocation(a_face,a_dx,probLo);
      const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
      Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
      Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
      CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);

// start hardwire
      Real wallThickness = 0.075;
      Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
      ParmParse pp;
      pp.get("wall_thickness",wallThickness);
      tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

      if (radius <= tubeRadius)
        {
          isInsideTube = true;
        }
    }
  bool isOutflow = (a_side==Side::Hi) && (a_idir==m_flowDir) && (!(isInsideTube));
 
//  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);

  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceGradPhi(a_faceFlux, a_face, a_comp,a_phi, a_probLo,
                               a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceGradPhi(a_faceFlux, a_face, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
}

//velocity stencil
void InflowOutflowHelmholtzDomainBC::getFluxStencil(      VoFStencil&      a_stencil,
                                                    const VolIndex&        a_vof,
                                                    const int&             a_comp,
                                                    const RealVect&        a_dx,
                                                    const int&             a_idir,
                                                    const Side::LoHiSide&  a_side,
                                                    const EBISBox&         a_ebisBox)
{
  bool isInsideTube = false;
  if ((m_doJet2) && (a_side==Side::Hi) && (a_idir == m_flowDir))
    {
      RealVect probLo = RealVect::Zero;
      RealVect loc = EBArith::getVofLocation(a_vof,a_dx,probLo);
      loc[a_idir] += 0.5*a_dx[a_idir]; // loc at face center
      const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
      Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
      Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
      CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);

// start hardwire
      Real wallThickness = 0.075;
      Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
      ParmParse pp;
      pp.get("wall_thickness",wallThickness);
      tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

      if (radius <= tubeRadius)
        {
          isInsideTube = true;
        }
    }
  bool isOutflow = (a_side==Side::Hi) && (a_idir==m_flowDir) && (!(isInsideTube)); 

  ParmParse pp;
  int bcType; // 0 for wall and 1 for free stream
  pp.get("bcType", bcType);

  if (bcType > 1) MayDay::Error("wrong BC type");

//  bool isOutflow = (a_side==Side::Hi) && (a_idir==m_flowDir);
  bool isSlipWall= ((bcType == 0) && (a_idir!=m_flowDir) && (a_idir != a_comp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));

  bool isFreeStreamBndry = ((bcType == 1) && (a_idir != m_flowDir));

  if (isOutflow || isSlipWall || isFreeStreamBndry)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setEBOrder(2);
      diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
}

///
/**
   This never gets called.  InflowOutflowPoissonDomainBC::getFaceVel takes care of it.
*/
void
InflowOutflowHelmholtzDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  MayDay::Error("InflowOutflowHelmholtzDomainBC::getFaceVel is not needed for Helmholtz.");
}
///
/**
   This is called by EBAMRPoissonOp::applyDomainFlux in EBAMRPoissonOp::applyOp
   for reg cells.
   For boundary conditions on velocity in viscous operator.
*/
void InflowOutflowHelmholtzDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                                 const BaseFab<Real>&  a_phi,
                                                 const RealVect&       a_probLo,
                                                 const RealVect&       a_dx,
                                                 const int&            a_idir,
                                                 const Side::LoHiSide& a_side,
                                                 const DataIndex&      a_dit,
                                                 const Real&           a_time,
                                                 const bool&           a_useHomogeneous)
{
  //vel: outflow is neumann. all others dirichlet.  inflow uses inflow vel as the value
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  
  ParmParse pp;
  int bcType;
  pp.get("bcType", bcType); // 0 solid wall; 1 free stream
  if (bcType > 1) MayDay::Error("wrong BC type");

  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((bcType == 0) && (a_idir!=m_flowDir) && (a_idir != velcomp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));
  bool isFreeStreamBndry = ((bcType == 1) && (a_idir != m_flowDir));
//  bool isVelNeum =  (isOutflow || isSlipWall);

  if (isSlipWall || isFreeStreamBndry)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
  else if (isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if (velcomp==m_flowDir)
        {
          if (!m_doJet1PoiseInflow)
            {
              diriBC.setValue(m_jet1inflowVel);
            }
          else if (m_doJet1PoiseInflow)
            {
              diriBC.setFunction(m_jet1PoiseInflowFunc);
            }
        }
      else
        {
          diriBC.setValue(0.0);
        }

      //basefab flux--EBAMRPoissonOp::applyDomainFlux calls this directly for viscous operator
      if (s_higherOrderHelmBC)
        {
          diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else
        {
          diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
    }
  else if (isOutflow)
    {
      if (!m_doJet2)
        {
          NeumannPoissonDomainBC neumannBC;
          neumannBC.setValue(0.);
          neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else if (m_doJet2)
        {
          const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
          Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
          CH_assert(a_phi.nComp() == 1);

// start hardwire
          Real wallThickness = 0.075;
          Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
          ParmParse pp;
          pp.get("wall_thickness",wallThickness);
          tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

          for (int comp=0; comp<a_phi.nComp(); comp++)
            {
              const Box& box = a_faceFlux.box();
              int iside = -1; //a_side == Side::Hi
              BoxIterator bit(box);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  RealVect loc = EBArith::getIVLocation(iv,a_dx,a_probLo);
                  loc[a_idir] -= iside * 0.5 * a_dx[a_idir];//point is now at the face center     
                  CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);
                  bool isInsideTube = false;
                  Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
                  if (radius <= tubeRadius)
                    {
                      isInsideTube = true;
                    }
                  if (isInsideTube) // Dirichlet
                    {
                      Real value = 0.0; // takes care of tangential comps
                      if (velcomp == m_flowDir)
                        {
                          if (m_doJet2PoiseInflow)
                            {
                              value = m_jet2PoiseInflowFunc->getVel(radius)[m_flowDir];
                            }
                          else if (!(m_doJet2PoiseInflow))
                            {
                              value = m_jet2inflowVel;
                            }
                        }
                      if (s_higherOrderHelmBC)
                        {
                          IntVect iv1 = bit();
                          Real phi0 = a_phi(iv1,comp);
                          iv1[a_idir] += iside;
                          Real phi1 = a_phi(iv1,comp);
                          iv1[a_idir] -= iside;
                          a_faceFlux(iv,comp) = -iside*(8.0*value + phi1 - 9.0*phi0)/(3.0*a_dx[a_idir]);
                        }
                      else if (!s_higherOrderHelmBC)
                        {
                          Real ihdx = 2.0 / a_dx[a_idir];
                          Real phiVal = a_phi(iv,comp);
                          a_faceFlux(iv,comp) = iside * ihdx * (phiVal - value);
                        }
                    }
                  else if (!isInsideTube) // Neumann
                    {
                      Real value = 0.0;
                      a_faceFlux(iv,comp) = iside * value;
                    }
                }
            } 
        }
    }
  else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
        {
          diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else
        {
          diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
    }
}

//called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
void InflowOutflowHelmholtzDomainBC::getFaceFlux(Real&                 a_faceFlux,
                                                 const VolIndex&       a_vof,
                                                 const int&            a_comp,
                                                 const EBCellFAB&      a_phi,
                                                 const RealVect&       a_probLo,
                                                 const RealVect&       a_dx,
                                                 const int&            a_idir,
                                                 const Side::LoHiSide& a_side,
                                                 const DataIndex&      a_dit,
                                                 const Real&           a_time,
                                                 const bool&           a_useHomogeneous)
{
  //vel: outflow is Neumann. all others Dirichlet
  int velcomp =  DirichletPoissonEBBC::s_velComp;

  ParmParse pp;
  int bcType;
  pp.get("bcType", bcType);
  if (bcType > 1) MayDay::Error("wrong BC type");

  bool isFreeStreamBndry = ((bcType == 1) && (a_idir != m_flowDir));
  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((bcType == 0) && (a_idir!=m_flowDir) && (a_idir != velcomp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));
//  bool isVelNeum =  (isOutflow || isSlipWall);

  if (isSlipWall || isFreeStreamBndry)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if (isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if (velcomp==m_flowDir)
        {
          if (!m_doJet1PoiseInflow)
          {
            diriBC.setValue(m_jet1inflowVel);
          }
          else if (m_doJet1PoiseInflow)
            {
              diriBC.setFunction(m_jet1PoiseInflowFunc);
            }
        }
      else
        {
          diriBC.setValue(0.0);
        }
      //called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
      if (s_higherOrderHelmBC)
        {
          diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else
        {
          diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        }
    }
  else if (isOutflow)
    {
      if (!m_doJet2)
        {
          NeumannPoissonDomainBC neumannBC;
          neumannBC.setValue(0.);
          neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        }
      else if (m_doJet2)
        {
          const IntVect& iv = a_vof.gridIndex();
          RealVect loc = EBArith::getIVLocation(iv,a_dx,a_probLo);
          loc[a_idir] += 0.5*a_dx[a_idir]; // loc is face centered
          const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
          Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
          CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);
          bool isInsideTube = false;
          Real radius = m_jet2PoiseInflowFunc->getRadius(loc);

// start hardwire
          Real wallThickness = 0.075;
          Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
          ParmParse pp;
          pp.get("wall_thickness",wallThickness);
          tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

          if (radius <= tubeRadius)
            {
              isInsideTube = true;
            }

          if (!(isInsideTube))
            {
              NeumannPoissonDomainBC neumannBC;
              neumannBC.setValue(0.);
              neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                    a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
            }
          else if (isInsideTube)
            {
              DirichletPoissonDomainBC diriBC;
              if (velcomp==m_flowDir)
                {
                  if (!m_doJet2PoiseInflow)
                    {
                      diriBC.setValue(m_jet2inflowVel);
                    }
                  else if (m_doJet2PoiseInflow)
                    {
                      diriBC.setFunction(m_jet2PoiseInflowFunc);
                    }
                }
              else
                {
                  diriBC.setValue(0.0);
                }
              //called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
              if (s_higherOrderHelmBC)
                {
                  diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
                }
              else
                {
                  diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
                }
            }
        }
    }
  else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
        {
          diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else
        {
          diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        }
    }
}

//////called by Helm solver for domain bc on phi
bool InflowOutflowHelmholtzDomainBC::
isDirichletDom(const VolIndex&   a_ivof,
               const VolIndex&   a_jvof,
               const EBCellFAB&  a_phi) const
{
/*
  FaceIndex face(a_ivof,a_jvof);
  int side = face.faceSign(a_ivof);
  int idir = face.direction();
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  bool isOutflow =  ((side== 1) && (idir==m_flowDir));
  bool isSlipWall = ((idir!=m_flowDir) && (idir != velcomp) && ((m_doSlipWallsHi[idir]==1 && side == 1)||(m_doSlipWallsLo[idir]==1 && side == -1)));
  bool isVelNeum =  (isOutflow || isSlipWall);
  if (isVelNeum)
    {
      return false;
    }
  else
    {
      return true;
    }
*/
  MayDay::Error("InflowOutflowPoisonDomainBC::isDirichletDom not set");
  return true; // to shut the compiler
}

//called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
void InflowOutflowHelmholtzDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
                                                      const VolIndex&       a_vof,
                                                      const int&            a_comp,
                                                      const EBCellFAB&      a_phi,
                                                      const RealVect&       a_probLo,
                                                      const RealVect&       a_dx,
                                                      const int&            a_idir,
                                                      const Side::LoHiSide& a_side,
                                                      const DataIndex&      a_dit,
                                                      const Real&           a_time)
{
  //vel: outflow is Neumann. all others Dirichlet
  int velcomp =  DirichletPoissonEBBC::s_velComp;

  ParmParse pp;
  int bcType;
  pp.get("bcType", bcType);
  if (bcType > 1) MayDay::Error("wrong BC type");

  bool isFreeStreamBndry = ((bcType == 1) && (a_idir!=m_flowDir));
  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((bcType == 0) && (a_idir!=m_flowDir) && (a_idir != velcomp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));
//  bool isVelNeum =  (isOutflow || isSlipWall);

  if (isSlipWall || isFreeStreamBndry)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit,a_time);
    }
  else if (isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if (velcomp==m_flowDir)
        {
          if (!m_doJet1PoiseInflow)
          {
            diriBC.setValue(m_jet1inflowVel);
          }
          else if (m_doJet1PoiseInflow)
            {
              diriBC.setFunction(m_jet1PoiseInflowFunc);
            }
        }
      else
        {
          diriBC.setValue(0.0);
        }
      //called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
      if (s_higherOrderHelmBC)
        {
          diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time);
        }
      else
        {
          diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time);
        }
    }
  else if (isOutflow)
    {
      if (!m_doJet2)
        {
          NeumannPoissonDomainBC neumannBC;
          neumannBC.setValue(0.);
          neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                     a_dx, a_idir, a_side, a_dit,a_time);
        }
      else if (m_doJet2)
        {
          const IntVect& iv = a_vof.gridIndex();
          RealVect loc = EBArith::getIVLocation(iv,a_dx,a_probLo);
          loc[a_idir] += 0.5*a_dx[a_idir]; // loc is face centered
          const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
          Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
          CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);
          bool isInsideTube = false;
          Real radius = m_jet2PoiseInflowFunc->getRadius(loc);

// start hardwire
          Real wallThickness = 0.075;
          Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
          ParmParse pp;
          pp.get("wall_thickness",wallThickness);
          tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire

          if (radius <= tubeRadius)
            {
              isInsideTube = true;
            }

          if (!(isInsideTube))
            {
              NeumannPoissonDomainBC neumannBC;
              neumannBC.setValue(0.);
              neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                         a_dx, a_idir, a_side, a_dit,a_time);
            }
          else if (isInsideTube)
            {
              DirichletPoissonDomainBC diriBC;
              if (velcomp==m_flowDir)
                {
                  if (!m_doJet2PoiseInflow)
                    {
                      diriBC.setValue(m_jet2inflowVel);
                    }
                  else if (m_doJet2PoiseInflow)
                    {
                      diriBC.setFunction(m_jet2PoiseInflowFunc);
                    }
                }
              else
                {
                  diriBC.setValue(0.0);
                }
              //called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
              if (s_higherOrderHelmBC)
                {
                  diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time);
                }
              else
                {
                  diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time);
                }
            } 
        }
    }
  else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
        {
          diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time);
        }
      else
        {
          diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time);
        }
    }
}

void InflowOutflowHelmholtzDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                    const FaceIndex&      a_face,
                                                    const int&            a_comp,
                                                    const EBCellFAB&      a_phi,
                                                    const RealVect&       a_probLo,
                                                    const RealVect&       a_dx,
                                                    const int&            a_idir,
                                                    const Side::LoHiSide& a_side,
                                                    const DataIndex&      a_dit,
                                                    const Real&           a_time,
                                                    const bool&           a_useAreaFrac,
                                                    const RealVect&       a_centroid,
                                                    const bool&           a_useHomogeneous)
{
  MayDay::Error("InflowOutflowHelmholtzDomainBC::getFaceGradPhi not needed");
}

#include "NamespaceFooter.H"
