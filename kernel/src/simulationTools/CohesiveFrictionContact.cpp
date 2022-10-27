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

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

using namespace RELATION;

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb, int numericsSolverId):
  FrictionContact(dimPb, SP::SolverOptions(solver_options_create(numericsSolverId),
                  solver_options_delete))
{
  if(! _q_cohesion)
    _q_cohesion.reset(new SiconosVector(LinearOSNS::maxSize()));

  if (_assemblyType == REDUCED_BLOCK or  _assemblyType == REDUCED_DIRECT)
  {
    if(! _V)
    {
      switch(_numericsMatrixStorageType)
      {
      case NM_DENSE:
      case NM_SPARSE:
      {
        _V.reset(new OSNSMatrix(0, _numericsMatrixStorageType));
        break;
      }
      case NM_SPARSE_BLOCK:
      {
        // = number of Interactionin the largest considered indexSet
        if(indexSetLevel() != LEVELMAX && simulation()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() > indexSetLevel())
        {
          _V.reset(new OSNSMatrix(simulation()->indexSet(indexSetLevel())->size(), _numericsMatrixStorageType));
        }
        else
        {
          _V.reset(new OSNSMatrix(1, _numericsMatrixStorageType));
        }
        break;
      }
      {
        default:
          THROW_EXCEPTION("LinearOSNS::initOSNSMatrix unknown _storageType");
      }
      }
    }
  }
}

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb, SP::SolverOptions options):
  FrictionContact(dimPb, options)
{
  if(! _q_cohesion)
    _q_cohesion.reset(new SiconosVector(LinearOSNS::maxSize()));

  if (_assemblyType == REDUCED_BLOCK or  _assemblyType == REDUCED_DIRECT)
  {
    if(! _V)
    {
      switch(_numericsMatrixStorageType)
      {
      case NM_DENSE:
      case NM_SPARSE:
      {
        _V.reset(new OSNSMatrix(0, _numericsMatrixStorageType));
        break;
      }
      case NM_SPARSE_BLOCK:
      {
        // = number of Interactionin the largest considered indexSet
        if(indexSetLevel() != LEVELMAX && simulation()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() > indexSetLevel())
        {
          _V.reset(new OSNSMatrix(simulation()->indexSet(indexSetLevel())->size(), _numericsMatrixStorageType));
        }
        else
        {
          _V.reset(new OSNSMatrix(1, _numericsMatrixStorageType));
        }
        break;
      }
      {
        default:
          THROW_EXCEPTION("LinearOSNS::initOSNSMatrix unknown _storageType");
      }
      }
    }
  }



}


void CohesiveFrictionContact::compute_q_cohesion_Block(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)
{
  DEBUG_BEGIN("CohesiveFrictionContact::compute_q_cohesion_Block(SP::Interaction inter, unsigned int pos)\n");
  SP::InteractionsGraph indexSet = simulation()->indexSet(0);

  OneStepIntegrator& osi1 = *indexSet->properties(vertex_inter).osi1;
  OneStepIntegrator& osi2 = *indexSet->properties(vertex_inter).osi2;

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
      // This has already been done for indexSet 1
      // we redo it for indexSet0. We should be more efficient
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs_cohesion = *(*indexSet->properties(vertex_inter).workVectors)[MoreauJeanOSI::OSNSP_RHS_COHESION];
      setBlock(osnsp_rhs_cohesion, _q_cohesion, sizeY, 0, pos);
    }
    else
      THROW_EXCEPTION("CohesiveFrictionContact::compute_q_cohesion_Block not yet implemented for OSI1 and OSI2 of type " + std::to_string(osi1Type)  + std::to_string(osi2Type));
  }
  DEBUG_EXPR(_q_cohesion->display());
  DEBUG_END("CohesiveFrictionContact::compute_q_cohesion_Block (SP::Interaction inter, unsigned int pos)\n");
}



void CohesiveFrictionContact::computeq(double time)
{
  LinearOSNS::computeq(time);

  // === Get index set from Simulation ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::InteractionsGraph indexSet0 = simulation()->indexSet(0);
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  if(_q_cohesion->size() != _sizeOutput_cohesion)
    _q_cohesion->resize(_sizeOutput_cohesion);
  _q_cohesion->zero();

  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    pos = indexSet0->properties(*ui).absolute_position;
    SP::Interaction inter = indexSet0->bundle(*ui);
    compute_q_cohesion_Block(*ui, pos);
  }


  DEBUG_EXPR(_q_cohesion->display(););


  //DEBUG_EXPR(_M->display(););
  NumericsMatrix * NM  = &*(_V->numericsMatrix());

  DEBUG_PRINT("before");
  DEBUG_EXPR(_q->display(););
  NM_gemv(1.0,
          NM,
          &*_q_cohesion->getArray(),
          1.0,
          &*_q->getArray());

  // DEBUG_EXPR(
  //   SP::SiconosVector q_add(new SiconosVector(*_q_cohesion));
  //   q_add->zero();
  //   NM_gemv(1.0,
  //           NM,
  //           &*_q_cohesion->getArray(),
  //           1.0,
  //           &*q_add->getArray());
  //   q_add->display();
  //   );
  DEBUG_EXPR(_q->display());
}
void CohesiveFrictionContact::computeV()
{
  if (_assemblyType == REDUCED_BLOCK)
  {

    InteractionsGraph& indexSet0 = *simulation()->indexSet(0);
    InteractionsGraph& indexSet1 = *simulation()->indexSet(1);
    indexSet0.update_vertices_indices();
    indexSet0.update_edges_indices();
    // Computes new _interactionBlocks if required
    updateInteractionBlocks(indexSet0);

    // _M->fillM(indexSet0, !_hasBeenUpdated);
    //  DEBUG_EXPR( _M->display(););

    // _V->fillM(indexSet0, !_hasBeenUpdated);
    // DEBUG_PRINT("complete V");
    // DEBUG_EXPR( _V->display(););

    _V->fillV(indexSet1, indexSet0, !_hasBeenUpdated);
    DEBUG_PRINT("partial V");
    DEBUG_EXPR( _V->display(););


  }
  else
    THROW_EXCEPTION("CohesiveFrictionContact::computeV unknown _assemblyTYPE");


  DEBUG_EXPR(_M->display(););
  // NumericsMatrix *   M_NM = _M->numericsMatrix().get();
  // if (M_NM )
  //   NM_display(M_NM);

  // getchar();

}

bool CohesiveFrictionContact::preCompute(double time)
{
  DEBUG_BEGIN(" CohesiveFrictionContact::preCompute(double time)\n");

  // Now we compute _V
  computeV();
  _sizeOutput_cohesion = _V->sizeColumn();
  DEBUG_PRINTF("_sizeOutput_cohesion = %i \n", _sizeOutput_cohesion );

  // _M and _q are computed on indexSet 1
  LinearOSNS::preCompute(time);

  unsigned int sizeInputIndexSet0 = simulation()->indexSet(0)->size();
  unsigned int sizeInputIndexSet1 = simulation()->indexSet(1)->size();
  DEBUG_PRINTF("sizeInputIndexSet0 = %i\t, sizeInputIndexSet1 = %i\t, _sizeOutput = %i\n", sizeInputIndexSet0, sizeInputIndexSet1, _sizeOutput);


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
  DEBUG_END(" CohesiveFrictionContact::preCompute(double time)\n");
  return true;
}

void CohesiveFrictionContact::postCompute()
{
  DEBUG_BEGIN("void CohesiveFrictionContact::postCompute()\n");
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  //DEBUG_EXPR(display());

  //*_z = *_z + *_q_cohesion;

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
bool CohesiveFrictionContact::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{

  float type_number= (float) (Type::value(nslaw) + 0.1 * nslaw.size());
  _nslawtype.insert(type_number);

  if (not (Type::value(nslaw) == Type::CohesiveZoneModelNIFNSL ||
	   Type::value(nslaw) == Type::NewtonImpactFrictionNSL  ))

  {
    THROW_EXCEPTION("\nCohesiveFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with CohesiveFrictionalContact one step nonsmooth problem. \n\
                      Compatible NonSmoothLaw are: CohesiveZoneModelNIFNSL (2D or 3D) \n");
    return false;
  }

  return true;
}
