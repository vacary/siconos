/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "CohesiveFrictionContact.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "CohesiveZoneModelNIFNSL.hpp"
#include "OSNSMatrix.hpp"
#include "NumericsMatrix.h"

#define DEBUG_NOCOLOR
#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

using namespace RELATION;

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb, int numericsSolverId):
  FrictionContact(dimPb, SP::SolverOptions(solver_options_create(numericsSolverId),
                  solver_options_delete))
{
  if(! _q_cohesion)
    _q_cohesion.reset(new SiconosVector(LinearOSNS::maxSize()));
}

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb, SP::SolverOptions options):
  FrictionContact(dimPb, options)
{
  if(! _q_cohesion)
    _q_cohesion.reset(new SiconosVector(LinearOSNS::maxSize()));
}


void CohesiveFrictionContact::compute_q_cohesion_Block(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)
{
  DEBUG_BEGIN("CohesiveFrictionContact::compute_q_cohesion_Block(SP::Interaction inter, unsigned int pos)\n");
  SP::InteractionsGraph indexSet = simulation()->indexSet(0);

  // // At most 2 DS are linked by an Interaction
  SP::DynamicalSystem ds1;
  SP::DynamicalSystem ds2;
  // --- Get the dynamical system(s) (edge(s)) connected to the current interaction (vertex) ---
  if(indexSet->properties(vertex_inter).source != indexSet->properties(vertex_inter).target)
  {
    DEBUG_PRINT("a two DS Interaction\n");
    ds1 = indexSet->properties(vertex_inter).source;
    ds2 = indexSet->properties(vertex_inter).target;
  }
  else
  {
    DEBUG_PRINT("a single DS Interaction\n");
    ds1 = indexSet->properties(vertex_inter).source;
    ds2 = ds1;
    // \warning this looks like some debug code, but it gets executed even with NDEBUG.
    // may be compiler does something smarter, but still it should be rewritten. --xhub
    InteractionsGraph::OEIterator oei, oeiend;
    for(std::tie(oei, oeiend) = indexSet->out_edges(vertex_inter);
        oei != oeiend; ++oei)
    {
      // note : at most 4 edges
      ds2 = indexSet->bundle(*oei);
      if(ds2 != ds1)
      {
        assert(false);
        break;
      }
    }
  }
  assert(ds1);
  assert(ds2);

  DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  OneStepIntegrator& osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
  OneStepIntegrator& osi2 = *DSG0.properties(DSG0.descriptor(ds2)).osi;

  OSI::TYPES osi1Type = osi1.getType();
  OSI::TYPES osi2Type = osi2.getType();

  SP::Interaction inter = indexSet->bundle(vertex_inter);


  SP::NonSmoothLaw nslaw = inter->nonSmoothLaw();
  unsigned int sizeY = nslaw->size();
  SP::CohesiveZoneModelNIFNSL nslaw_CohesiveZoneModelNIFNSL(std::dynamic_pointer_cast<CohesiveZoneModelNIFNSL>(nslaw));
  if (nslaw_CohesiveZoneModelNIFNSL)
  {
    if((osi1Type == OSI::MOREAUJEANOSI  && osi2Type == OSI::MOREAUJEANOSI)||
       (osi1Type == OSI::MOREAUDIRECTPROJECTIONOSI && osi2Type == OSI::MOREAUDIRECTPROJECTIONOSI))
    {
      //osi1.computeFreeOutput(vertex_inter, this); This has already been done
      SiconosVector& osnsp_rhs_cohesion = *(*indexSet->properties(vertex_inter).workVectors)[MoreauJeanOSI::OSNSP_RHS_COHESION];
      setBlock(osnsp_rhs_cohesion, _q_cohesion, sizeY, 0, pos);
    }
    else
      THROW_EXCEPTION("CohesiveFrictionContact::compute_q_cohesion_Block not yet implemented for OSI1 and OSI2 of type " + std::to_string(osi1Type)  + std::to_string(osi2Type));
  }
  //DEBUG_EXPR(_q_cohesion->display());
  DEBUG_END("CohesiveFrictionContact::compute_q_cohesion_Block (SP::Interaction inter, unsigned int pos)\n");
}



void CohesiveFrictionContact::computeq(double time)
{
  LinearOSNS::computeq(time);

  // === Get index set from Simulation ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  if(_q_cohesion->size() != _sizeOutput)
    _q_cohesion->resize(_sizeOutput);
  _q_cohesion->zero();

  for(std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    pos = indexSet->properties(*ui).absolute_position;
    SP::Interaction inter = indexSet->bundle(*ui);
    DEBUG_EXPR(inter->nonSmoothLaw()->display(););
    compute_q_cohesion_Block(*ui, pos);
  }

  //DEBUG_EXPR(_M->display(););
  NumericsMatrix * NM  = &*(_M->numericsMatrix());


  DEBUG_EXPR(_q_cohesion->display(););
  DEBUG_EXPR("before"; _q->display(););
  NM_gemv(1.0,
          NM,
          &*_q_cohesion->getArray(),
          1.0,
          &*_q->getArray());

  DEBUG_EXPR(
    SP::SiconosVector q_add(new SiconosVector(*_q_cohesion));
    q_add->zero();
    NM_gemv(1.0,
            NM,
            &*_q_cohesion->getArray(),
            1.0,
            &*q_add->getArray());
    q_add->display();
    );
  DEBUG_EXPR(_q->display());
}

bool CohesiveFrictionContact::preCompute(double time)
{
  LinearOSNS::preCompute(time);
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  if(_keepLambdaAndYState)
  {
    InteractionsGraph::VIterator ui, uiend;
    for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      Interaction& inter = *indexSet.bundle(*ui);
      SP::NonSmoothLaw nslaw = inter.nonSmoothLaw();
      SP::CohesiveZoneModelNIFNSL nslaw_CohesiveZoneModelNIFNSL(std::dynamic_pointer_cast<CohesiveZoneModelNIFNSL>(nslaw));
      if (nslaw_CohesiveZoneModelNIFNSL)
      {
        // Get the position of inter-interactionBlock in the vector w
        // or z
        unsigned int pos = indexSet.properties(*ui).absolute_position;
        SiconosVector& osnsp_rhs_cohesion = *(*indexSet.properties(*ui).workVectors)[MoreauJeanOSI::OSNSP_RHS_COHESION];
        
        for (int k =0; k < osnsp_rhs_cohesion.size(); k++)
        {
          (*_z)(pos+k) -= osnsp_rhs_cohesion(k);
        }
      }
    }
  }
  return true;
}

void CohesiveFrictionContact::postCompute()
{
  DEBUG_BEGIN("void CohesiveFrictionContact::postCompute()\n");
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  //DEBUG_EXPR(display());

  *_z = *_z + *_q_cohesion;

  DEBUG_EXPR(_w->display(););
  DEBUG_EXPR(_z->display(););
  // === Get index set from Topology ===
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===

  unsigned int pos = 0;

  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet.bundle(*ui);
    // Get the  position of inter-interactionBlock in the vector w
    // or z
    pos = indexSet.properties(*ui).absolute_position;

    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy _w/_z values, starting from index pos into y/lambda.

    //setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved in y !!



    setBlock(*_z, lambda, lambda->size(), pos, 0);
    DEBUG_EXPR(lambda->display(););
  }

  DEBUG_END("void CohesiveFrictionContact::postCompute()\n");
}
