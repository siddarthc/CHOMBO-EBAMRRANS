#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBArith.H"
#include "InflowOutflowAdvectBC.H"
#include "DirichletPoissonEBBC.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

extern Real g_simulationTime;

/****************************/
void
InflowOutflowAdvectBC::
fluxBC(EBFluxFAB&            a_primGdnv,
       const EBCellFAB&      a_primCenter,
       const EBCellFAB&      a_primExtrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_dir)
{
  CH_assert(m_isDefined);
  CH_assert(!m_domain.isPeriodic(a_dir));

  //gets face centered region of data
  Box FBox = a_primGdnv[a_dir].getRegion();
  Box cellBox = a_faceBox;
  CH_assert(a_primGdnv[a_dir].nComp() == 1);

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  //figure out if we are on a wall and if its an inflow and the inflow value
  bool isInflow   = ((a_side == Side::Lo) && (a_dir==m_flowDir));
  bool isOutflow  = ((a_side == Side::Hi) && (a_dir==m_flowDir));

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box bndryBox = adjCellBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      bndryBox.shift(a_dir,-isign);

      IntVectSet ivs(bndryBox);
      for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
          for (int iface= 0; iface < bndryFaces.size(); iface++)
            {
              //this all will work for the case of spatially varying inflow
              //for inhomogeneous other faces, this is wrong.
              const FaceIndex& face = bndryFaces[iface];
              Real velComp = 1.e99;
              if (isInflow && (m_velComp == m_flowDir))
                {
                  if (!m_doJet1PoiseInflow)
                    {
                      velComp = m_jet1inflowVel;
                    }
                  else if (m_doJet1PoiseInflow)
                    {
                      //set coordinates based on slip walls
                      RealVect prob_lo = RealVect::Zero;
                      //assumes inflow boundary is always on lo side of domain
                      const RealVect loc  = EBArith::getFaceLocation(face, m_dx, prob_lo);
                      Real radius = m_jet1PoiseInflowFunc->getRadius(loc);
                      velComp = m_jet1PoiseInflowFunc->getVel(radius)[m_flowDir];
                    }
                }
              else if (isInflow && (m_velComp != m_flowDir))
                {
                  velComp = 0.0;
                }
              else if (isOutflow)
                {
                   if (!m_doJet2)
                     {
                       velComp = a_primExtrap(vof, 0);
                     }
                   else if (m_doJet2)
                     {
                       RealVect prob_lo = RealVect::Zero;
                       const RealVect tubeCenter = m_jet2PoiseInflowFunc->getTubeCenter();
                       Real tubeRadius = m_jet2PoiseInflowFunc->getTubeRadius();
                       bool isInsideTube = false;
                       const RealVect loc  = EBArith::getFaceLocation(face, m_dx, prob_lo);                        
                       CH_assert(loc[m_flowDir] == tubeCenter[m_flowDir]);

                       Real radius = m_jet2PoiseInflowFunc->getRadius(loc);
// start hardwire
                       Real wallThickness = 0.075;
                       ParmParse pp;
                       pp.get("wall_thickness",wallThickness);
                       Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
                       tubeRadius += wallThickness * bcFactor /2.0;
// end hardwire
                       if (radius <= tubeRadius)
                         {
                           isInsideTube = true;
                         }
                       if (isInsideTube && (m_velComp == m_flowDir))
                         {
                           if (m_doJet2PoiseInflow)
                             {
                               velComp = m_jet2PoiseInflowFunc->getVel(radius)[m_flowDir];
                             }
                           else if (!m_doJet2PoiseInflow)
                             {
                               velComp = m_jet2inflowVel;
                             }
                         }
                       else if (isInsideTube && (m_velComp != m_flowDir))
                         {
                           velComp = 0.0;
                         }
                       else if (!(isInsideTube))
                         {
                           velComp = a_primExtrap(vof, 0);
                         } 
                     }
                }
              else //solid wall or free stream 
                {
                  ParmParse pp;
                  int bcType; // 0 for solid wall, 1 for free stream
                  pp.get("bcType", bcType);

                  if (a_dir == m_velComp)
                    {
                      if (bcType == 0) // solid wall
                        {
                          velComp = 0.0;
                        }
                      else if (bcType == 1) // free stream
                        {
                          velComp = a_primExtrap(vof, 0); 
                        }
                      else 
                        {
                          MayDay::Error("wrong bcType");
                        }
                    }
                  else
                    {
                      if (bcType == 0) // solid wall
                        {  
                          velComp = a_primExtrap(vof, 0);
                        }
                      else if (bcType == 1) // free stream
                        {
                          velComp = a_primExtrap(vof, 0); 
                        }
                      else
                        {
                          MayDay::Error("wrong bcType"); 
                        }
                    }
                }
              a_primGdnv[a_dir](face, 0) = velComp;
            }
        }
    }
}

#include "NamespaceFooter.H"
