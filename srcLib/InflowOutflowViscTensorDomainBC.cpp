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
//  int velcomp =  DirichletViscousTensorEBBC::s_velComp;

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
          
        }
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

}

#include "NamespaceFooter.H"
