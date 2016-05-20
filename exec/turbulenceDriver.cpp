#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"


#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBFABView.H"
#include "memtrack.H"
#include "AMRINSUtils.H"
#include "CH_Attach.H"
#include "EBAMRNoSubcycle.H"
#include "InflowOutflowIBC.H"
#include "EBScalarAdvectBCFactory.H"
#include "PoiseuilleInflowBCValue.H"
#include "EBAMRTransportFactory.H"
#include "EBAMRTransportParams.H"
#include "EBPatchTransport.H"
#include "EBFABView.H"
#include <iostream>
#include "UsingNamespace.H"
#include "memusage.H"
#include "memtrack.H"

/***************/
void ebamriTransport(const AMRParameters& a_params,
                     const ProblemDomain& a_coarsestDomain)
{

  // begin define INS solver
  CH_TIMERS("ebamriTransport_driver");
  CH_TIMER("define_ebamrnosubcycle_solver", t3);
  CH_TIMER("init_ebamrnosubcycle_solver",   t4);
  CH_TIMER("run ebamrnosubcycle_solver",   t5);

  // read inputs
  ParmParse pp;

  int flowDir;
  pp.get("flow_dir", flowDir);
  Vector<int> nCells;
  pp.getarr("n_cell",  nCells,0,SpaceDim);

  Real jet1inflowVel;
  pp.get("jet1_inflow_vel", jet1inflowVel);

  int idoJet2;
  pp.get("do_jet2", idoJet2);
  bool doJet2 = idoJet2==1;

  Real jet2inflowVel;
  pp.get("jet2_inflow_vel", jet2inflowVel);

  Real viscosity = 0.0;
  pp.get("viscosity", viscosity);

  int idoSlipWalls;
  pp.get("do_slip_walls", idoSlipWalls);
  bool doSlip = idoSlipWalls==1;
  IntVect doSlipWallsLo = idoSlipWalls*IntVect::Unit;
  IntVect doSlipWallsHi = idoSlipWalls*IntVect::Unit;
  Vector<int> slipWallsLo,slipWallsHi;
  if(doSlip)
    {
      pp.getarr("do_slip_walls_hi",slipWallsHi,0,SpaceDim);
      pp.getarr("do_slip_walls_lo",slipWallsLo,0,SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          doSlipWallsHi[idir] = slipWallsHi[idir];
          doSlipWallsLo[idir] = slipWallsLo[idir];
        }
    }

  int orderEBBC = 1;
  pp.query("order_ebbc", orderEBBC);

  bool doJet1PoiseInflow = false;
  pp.query("jet1_poiseuille_inflow", doJet1PoiseInflow);

  bool doJet2PoiseInflow = false;
  pp.query("jet2_poiseuille_inflow", doJet2PoiseInflow);

  bool initPoiseData = false;
  pp.query("poiseuille_init", initPoiseData);
  if(initPoiseData)  
  pout() << "Doing Poiseuille initialization" << endl;

  RefCountedPtr<PoiseuilleInflowBCValue> jet1PoiseBCValue;//make this BaseBCValue if also doing constant values
  RefCountedPtr<PoiseuilleInflowBCValue> jet2PoiseBCValue;

//  if(doPoiseInflow)
//    {
//      pout() << "Doing Poiseuille inflow" << endl;
      RealVect jet1centerPt, tubeAxis;
      Real jet1tubeRadius;
      Vector<Real> jet1centerPtVect, tubeAxisVect;
      pp.get("jet1_poise_profile_radius", jet1tubeRadius);
      pp.getarr("jet1_poise_profile_center_pt",  jet1centerPtVect,0,SpaceDim);
      pp.getarr("poise_profile_axis",       tubeAxisVect,0,SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          jet1centerPt[idir] = jet1centerPtVect[idir];
          tubeAxis[idir] = tubeAxisVect[idir];
        }

      Real maxVelFactor;//= 1.5 for planar geometry, = 2.0 for cylindrical
      pp.get("poise_maxvel_factor", maxVelFactor);
      Real jet1maxVel = maxVelFactor*jet1inflowVel;

      jet1PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue>(new PoiseuilleInflowBCValue(jet1centerPt,tubeAxis,jet1tubeRadius,jet1inflowVel,flowDir,viscosity));
//    }

      RealVect jet2centerPt, jet2TubeOrig, jet2TubeEnd;
      Real jet2tubeRadius;
      Vector<Real> jet2centerPtVect, jet2TubeOrigVect, jet2TubeEndVect;
      pp.get("jet2_poise_profile_radius", jet2tubeRadius);
      pp.getarr("jet2_poise_profile_center_pt",  jet2centerPtVect,0,SpaceDim);
      pp.getarr("jet2_tube_origin",  jet2TubeOrigVect,0,SpaceDim);
      pp.getarr("jet2_tube_end",  jet2TubeEndVect,0,SpaceDim); 
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          jet2centerPt[idir] = jet2centerPtVect[idir];
          jet2TubeOrig[idir] = jet2TubeOrigVect[idir];
          jet2TubeEnd[idir] = jet2TubeEndVect[idir];
        }

      Real jet2maxVel = maxVelFactor*jet2inflowVel;

      jet2PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue>(new PoiseuilleInflowBCValue(jet2centerPt,tubeAxis,jet2tubeRadius,jet2inflowVel,flowDir,viscosity));

  InflowOutflowIBCFactory ibc(flowDir,
                              doJet2,
                              jet1inflowVel,
                              jet2inflowVel,
                              orderEBBC,
                              jet2TubeOrig,
                              jet2TubeEnd,
                              doSlipWallsHi,
                              doSlipWallsLo,
                              doJet1PoiseInflow,
                              doJet2PoiseInflow,
                              initPoiseData,
                              jet1PoiseBCValue,
                              jet2PoiseBCValue);

  CH_START(t3);

  EBAMRNoSubcycle kahuna(a_params, ibc, a_coarsestDomain, viscosity);

  CH_STOP(t3);
  
  // end define INS solver

  // begin define transport solver
  
  // 1. setup the ibc  
  Real scalarInflowValue = 0.0;
  bool doScalarInit = false;
  int iScalInit;
  pp.get("scalar_inflow_value", scalarInflowValue);
  pp.get("scal_init", iScalInit);
  doScalarInit = (iScalInit==1);
  RefCountedPtr<EBScalarAdvectBCFactory> advectBCFact = 
    RefCountedPtr<EBScalarAdvectBCFactory>
    (new EBScalarAdvectBCFactory(flowDir, scalarInflowValue, jet2TubeOrig, jet2TubeEnd, doScalarInit, jet2PoiseBCValue));

  EBAMRTransportParams transParams;
  // fillTransParams(transParams);
  Real schmidt;
  Real diffCoeff;
  pp.get("schmidt_number", schmidt);
  diffCoeff = viscosity/schmidt;
  transParams.m_diffusionCoeff = diffCoeff;

  bool scalSlopeLimiting = true;
  int ilimit;
  pp.get("limit_scal_slope", ilimit);
  scalSlopeLimiting = (ilimit==1);

  Real maxScalVal = 1.0;
  pp.get("max_scal_value", maxScalVal);

  Real minScalVal = 0.0;
  pp.get("min_scal_value", minScalVal);

  bool setScalMaxMin = true;
  int ilimit2;
  pp.get("set_scal_maxmin", ilimit2);
  setScalMaxMin = (ilimit2==1);

  int ifourth, iflatten, iartvisc;
  pp.get("scal_use_fourth_order_slopes", ifourth);
  pp.get("scal_use_flatten", iflatten);
  pp.get("scal_use_art_visc", iartvisc);

  bool useFourthOrderSlopes = (ifourth  == 1);
  bool useFlattening        = (iflatten == 1);
  bool useArtificialVisc    = (iartvisc == 1);

  RefCountedPtr<EBPhysIBCFactory> ibcFact = static_cast<RefCountedPtr<EBPhysIBCFactory> >(advectBCFact);

  RefCountedPtr<EBPatchTransportFactory> patchTransport = 
   RefCountedPtr<EBPatchTransportFactory>
   (new EBPatchTransportFactory(ibcFact, scalSlopeLimiting, maxScalVal, minScalVal, setScalMaxMin, useFourthOrderSlopes, useFlattening, useArtificialVisc));

  RefCountedPtr<EBPatchGodunovFactory> patchFact = static_cast<RefCountedPtr<EBPatchGodunovFactory> >(patchTransport);
  // 2. setup transport solver
  EBAMRTransportFactory transSolverFact(transParams, patchFact, true);  
  // end define transport solver

  kahuna.defineExtraSolver(&transSolverFact);
  

  CH_START(t4);
  if (!pp.contains("restart_file"))
    {
      pout() << "starting fresh AMR run" << endl;
      kahuna.setupForAMRRun();
    }
  else
    {
      std::string restart_file;
      pp.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;
      kahuna.setupForRestart(restart_file);
    }
  CH_STOP(t4);

  int maxStep;
  pp.get("max_step", maxStep);

  Real stopTime = 0.0;
  pp.get("max_time",stopTime);

  CH_START(t5);

  Real fixedDt = 0.;
  pp.query("fixed_dt", fixedDt);
  if(fixedDt > 1.e-12)
    {
      kahuna.useFixedDt(fixedDt);
    }

  kahuna.run(stopTime, maxStep);

  CH_STOP(t5);

}
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif
  //Scoping trick
  {
    CH_TIMERS("uber_timers");
    CH_TIMER("define_geometry", t1);
    CH_TIMER("run", t2);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    //Check for an input file
    char* inFile = NULL;
    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }

    //Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    ProblemDomain coarsestDomain;
    RealVect coarsestDx;
    AMRParameters params;
    getAMRINSParameters(params, coarsestDomain);
    int numFilt;
    pp.get("num_filter_iterations", numFilt);
    params.m_numFilterIterations = numFilt;

    int gphiIterations;
    pp.get("num_gphi_iterations", gphiIterations);
    params.m_gphiIterations = gphiIterations;

    int initIterations;
    pp.get("num_init_iterations", initIterations);
    params.m_initIterations = initIterations;

    bool doRegridSmoothing;
    pp.get("do_regrid_smoothing", doRegridSmoothing);
    params.m_doRegridSmoothing = doRegridSmoothing;

    CH_START(t1);

    //define geometry
    AMRINSGeometry(params, coarsestDomain);

    CH_STOP(t1);

    CH_START(t2);
    ebamriTransport(params, coarsestDomain);
    CH_STOP(t2);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }//end scoping trick
#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
