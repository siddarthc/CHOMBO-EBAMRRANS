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

#include "InflowOutflowViscTensorDomainBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletViscousTensorEBBC.H"
#include "NeumannViscousTensorEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"
#include "PoiseuilleInflowBCValue.H"
#include "SlipWallViscousTensorDomainBC.H"
#include "PoiseuilleInflowFuncEval.H"
#include "FreeStreamFuncEval.H"

#include "NamespaceHeader.H"

void InflowOutflowViscTensorDomainBC::
whereAMI(bool a_isInsideTube, const RealVect& a_loc)
{
  CH_assert(a_loc[m_flowDir] == (m_jet2PoiseInflowFunc->getTubeCenter())[m_flowDir]);
  a_isInsideTube = false;
  Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
 
// start hardwire
  Real wallThickness = 0.075;
  Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
  ParmParse pp;
  pp.get("wall_thickness",wallThickness);
  tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire
 
  Real radius = m_jet2PoiseInflowFunc->getRadius(a_loc);
  if (radius <= tubeRadius)
   {
     a_isInsideTube = true;
   } 
}

void InflowOutflowViscTensorDomainBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
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

  ParmParse pp;
  int bcType;
  pp.get("bcType", bcType); // 0 solid wall; 1 free stream
  // need to fix the free stream condition!
  if (bcType > 1) MayDay::Error("wrong BC type");

  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((bcType == 0) && (a_idir!=m_flowDir) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));
  bool isFreeStreamBndry = ((bcType == 1) && (a_idir != m_flowDir));

  if (isSlipWall)
    {
      SlipWallViscousTensorDomainBC slipBC;
      slipBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      slipBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }  

  else if (isFreeStreamBndry)
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      RefCountedPtr<BaseBCFuncEval> funcEval(new FreeStreamFuncEval (m_jet1inflowVel, m_flowDir));
      diriBC.setFunction(funcEval);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }

  else if (isInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
 
      if (!m_doJet1PoiseInflow)
        {
          RefCountedPtr<BaseBCFuncEval> funcEval(new FreeStreamFuncEval (m_jet1inflowVel, m_flowDir));
          diriBC.setFunction(funcEval);
        }
      else
        {
          RefCountedPtr<BaseBCFuncEval> funcEval(new PoiseuilleInflowFuncEval (m_jet1PoiseInflowFunc));
          diriBC.setFunction(funcEval); 
        }
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }

  else if (isOutflow)
    {
      if (!m_doJet2)
        {
          NeumannViscousTensorDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          neumannBC.setValue(0.);
          neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
        }
      else if (m_doJet2)
        {
          // the painful part
          // figure out where I am
          // If I am in the tube, Dirichlet else Neumann

          // Life would be so much easier if the EBBC version works
          //   Its just sad now
          // flux is eta(grad B + grad B^T) + lambda* I div B   on the face
          CH_assert(a_phi.nComp() == SpaceDim); 
          const FArrayBox& phiFAB = (const FArrayBox& )a_phi; 
          const Box& faceBox = a_faceFlux.box();
          int iside = -1; //a_side == Side::Hi, its -1 for PoissonOp

           RefCountedPtr<BaseBCFuncEval> funcEval;

          if (!m_doJet2PoiseInflow)
            {
              funcEval = RefCountedPtr<BaseBCFuncEval>(new FreeStreamFuncEval (m_jet2inflowVel, m_flowDir));
            }
          else
            {
              funcEval = RefCountedPtr<BaseBCFuncEval>(new PoiseuilleInflowFuncEval (m_jet2PoiseInflowFunc));
            }

          int ncompPhi = SpaceDim;
          int ncompGph = SpaceDim*SpaceDim;
          FArrayBox faceDiv(faceBox, 1);
          FArrayBox faceGph(faceBox, ncompGph);
          FArrayBox phiValu(faceBox, ncompPhi);
          faceDiv.setVal(0.);
          faceGph.setVal(0.);

          BoxIterator bit(faceBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              RealVect loc = EBArith::getIVLocation(iv,a_dx,a_probLo);
              loc[a_idir] += iside * 0.5 * a_dx[a_idir];// difference between PiissonOp and this is that here the incoming box is half shifted
              bool isInsideTube;
              whereAMI(isInsideTube, loc);

              if (isInsideTube)
                {
                  for (int comp = 0; comp < SpaceDim; comp++)
                    {
                      if (!a_useHomogeneous)
                        {
                          Real value = funcEval->value(loc, comp, a_time);
                          phiValu(iv, comp) = value;
                          for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
                            {
                              if (derivDir != a_idir)
                                {
                                  Real deriv   = funcEval->derivative(loc, comp, derivDir);
                                  int gradComp = TensorCFInterp::gradIndex(comp, derivDir);
                                  faceGph(iv, gradComp) = deriv;
                                }  // end if
                            } // end derivDir
                        }  // end homogeneous
                      else
                        {
                          phiValu(iv, comp) = 0.;
                          faceGph(iv, comp) = 0.;
                        } // end homogeneous
                    } // end comp

                  // want to use the value for normal derivs
                  int divDir = a_idir;
                  // for now no higher order BC
                  Real ihdx = 2.0/a_dx[0];
                  for (int velDir = 0; velDir < SpaceDim; velDir++)
                    {
                      int gradcomp = TensorCFInterp::gradIndex(velDir,divDir);
                      faceGph(iv, gradcomp) = iside * ihdx * (phiFAB(iv, velDir) - phiValu(iv, velDir)); 
                    }

                } //end isInsideTube
              else // outflow
                {
                  // just homogeneous neumann 0
                  Real value = 0.;
                  for (int comp = 0; comp < SpaceDim; comp++)
                    {
                      faceGph(iv, comp) = 0.;
                    }
                } // end outflow 
            } // end boxIterator 
          getFluxFromGrad(a_faceFlux, faceGph, a_dit, a_idir); 
        } // end m_doJet2
    }
  else // no slip wall
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.setValue(0.);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);     
    }
}

void InflowOutflowViscTensorDomainBC::
getFaceFlux(Real&                 a_faceFlux,
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
  ParmParse pp;
  int bcType;
  pp.get("bcType", bcType); // 0 solid wall; 1 free stream
  // need to fix the free stream condition!
  if (bcType > 1) MayDay::Error("wrong BC type");

  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((bcType == 0) && (a_idir!=m_flowDir) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));
  bool isFreeStreamBndry = ((bcType == 1) && (a_idir != m_flowDir));

  if (isSlipWall)
    {
      SlipWallViscousTensorDomainBC slipbc;
      slipbc.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      slipbc.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }

  else if (isFreeStreamBndry)
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      RefCountedPtr<BaseBCFuncEval> funcEval(new FreeStreamFuncEval (m_jet1inflowVel, m_flowDir));
      diriBC.setFunction(funcEval);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }

  else if (isInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);

      if (!m_doJet1PoiseInflow)
        {
          RefCountedPtr<BaseBCFuncEval> funcEval(new FreeStreamFuncEval (m_jet1inflowVel, m_flowDir));
          diriBC.setFunction(funcEval);
        }
      else
        {
          RefCountedPtr<BaseBCFuncEval> funcEval(new PoiseuilleInflowFuncEval (m_jet1PoiseInflowFunc));
          diriBC.setFunction(funcEval);
        }

      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);    
    }

  else if (isOutflow)
    {
      if (!m_doJet2)
        {
          NeumannViscousTensorDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          neumannBC.setValue(0.0);
          neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                                a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        }
      else
        {
          const IntVect& iv = a_vof.gridIndex();
          RealVect loc = EBArith::getIVLocation(iv,a_dx,a_probLo);
          loc[a_idir] += 0.5*a_dx[a_idir]; // loc is face centered
          bool isInsideTube;
          whereAMI(isInsideTube, loc);
          if (isInsideTube)
            {
              RefCountedPtr<BaseBCFuncEval> funcEval;

              if (!m_doJet2PoiseInflow)
                {
                  funcEval = RefCountedPtr<BaseBCFuncEval>(new FreeStreamFuncEval (m_jet2inflowVel, m_flowDir));
                }
              else
                {
                  funcEval = RefCountedPtr<BaseBCFuncEval>(new PoiseuilleInflowFuncEval (m_jet2PoiseInflowFunc));
                }

              DirichletViscousTensorDomainBC diriBC;
              diriBC.setFunction(funcEval);
              diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
              diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
            }
          else
            {
              NeumannViscousTensorDomainBC neumannBC;
              neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
              neumannBC.setValue(0.0);
              neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                                    a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
            }
        }
    }

  else
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.setValue(0.);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

#include "NamespaceFooter.H"
