/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "TimeStepping.hpp"
#include "Topology.hpp"
#include "LCP.hpp"
#include "Relay.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
//#include "Interaction.hpp"
#include "OneStepIntegrator.hpp"
#include "EulerMoreauOSI.hpp"
#include "Interaction.hpp"
#include "EventsManager.hpp"
#include "FrictionContact.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "NonSmoothLaw.hpp"
#include "TypeName.hpp"
#include "Relation.hpp"
#include "BlockVector.hpp"
#include "NewtonEulerR.hpp"
#include "FirstOrderR.hpp"

#include <SiconosConfig.h>
#include <functional>
using namespace std::placeholders;

// #define DEBUG_BEGIN_END_ONLY
#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include <debug.h>

using namespace RELATION;

/** Pointer to function, used to set the behavior of simulation when
    ns solver failed.  If equal to null, use DefaultCheckSolverOutput
    else (set with setCheckSolverFunction) call the pointer below).
    Note FP: (temporary) bad method to set checkSolverOutput but it
    works ... It may be better to use plug-in?
*/
static CheckSolverFPtr checkSolverOutput = nullptr;

TimeStepping::TimeStepping(SP::NonSmoothDynamicalSystem nsds,
                           SP::TimeDiscretisation td,
                           SP::OneStepIntegrator osi,
                           SP::OneStepNSProblem osnspb)
  : Simulation(nsds,td), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonNbIterations(0),
    _newtonCumulativeNbIterations(0), _newtonOptions(SICONOS_TS_NONLINEAR),
    _newtonResiduDSMax(0.0), _newtonResiduYMax(0.0), _newtonResiduRMax(0.0),
    _computeResiduY(false),_computeResiduR(false),
    _isNewtonConverge(false),
    _newtonUpdateInteractionsPerIteration(false),_displayNewtonConvergence(false),
    _warnOnNonConvergence(true),
    _resetAllLambda(true)
{

  if(osi) insertIntegrator(osi);
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TS);
  if(osnspb) insertNonSmoothProblem(osnspb, SICONOS_OSNSP_TS_VELOCITY);
  SP::EulerMoreauOSI euosi(std::dynamic_pointer_cast<EulerMoreauOSI>(osi));
  if(euosi)
  {
    _computeResiduY = true;
    _computeResiduR = true;
  }

}

TimeStepping::TimeStepping(SP::NonSmoothDynamicalSystem nsds, SP::TimeDiscretisation td, int nb)
  : Simulation(nsds,td), _newtonTolerance(1e-6), _newtonMaxIteration(50), _newtonNbIterations(0),
    _newtonCumulativeNbIterations(0), _newtonOptions(SICONOS_TS_NONLINEAR),
    _newtonResiduDSMax(0.0), _newtonResiduYMax(0.0), _newtonResiduRMax(0.0), _computeResiduY(false),
    _computeResiduR(false),
    _isNewtonConverge(false),
    _newtonUpdateInteractionsPerIteration(false),_displayNewtonConvergence(false),
    _warnOnNonConvergence(true),
    _resetAllLambda(true)
{
  (*_allNSProblems).resize(nb);
}

// --- Destructor ---
TimeStepping::~TimeStepping()
{
}

void TimeStepping::insertIntegrator(SP::OneStepIntegrator osi)
{
  _allOSI->insert(osi);
  SP::EulerMoreauOSI euosi(std::dynamic_pointer_cast<EulerMoreauOSI>(osi));
  if(euosi)
  {
    _computeResiduY = true;
    _computeResiduR = true;
  }
}
// bool TimeStepping::predictorDeactivate(SP::Interaction inter, unsigned int i)
// {
//   double h = timeStep();
//   double y = inter->getYRef(i-1); // for i=1 it is the position -> historic notation y
//   double yDot = inter->getYRef(i); // for i=1 it is the velocity -> historic notation yDot
//   DEBUG_PRINTF("TS::predictorDeactivate yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y+0.5*h*yDot);
//   y += 0.5*h*yDot;
//   assert(!std::isnan(y));
//   if(y>0)
//     DEBUG_PRINTF("TS::predictorDeactivate DEACTIVATE.\n");
//   return (y>0);
// }

// bool TimeStepping::predictorActivate(SP::Interaction inter, unsigned int i)
// {
//   double h = timeStep();
//   double y = inter->getYRef(i-1); // for i=1 it is the position -> historic notation y
//   double yDot = inter->getYRef(i); // for i=1 it is the velocity -> historic notation yDot
//   DEBUG_PRINTF("TS::predictorActivate yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y+0.5*h*yDot);
//   y += 0.5*h*yDot;
//   assert(!std::isnan(y));
//   if(y<=0)
//     DEBUG_PRINTF("TS::predictorActivate ACTIVATE.\n");
//   return (y<=0);
// }

void TimeStepping::updateIndexSet(unsigned int i)
{
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.
  // boost::default_color_type is used to organize update in InteractionsGraph:
  // - white_color : undiscovered vertex (Interaction)
  // - gray_color : discovered vertex (Interaction) but searching descendants
  // - black_color : discovered vertex (Interaction) together with the descendants

  assert(_nsds);
  assert(_nsds->topology());

  SP::Topology topo = _nsds->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0 &&
         "TimeStepping::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  assert(indexSet0);
  assert(indexSet1);
  DynamicalSystemsGraph& DSG0= *nonSmoothDynamicalSystem()->dynamicalSystems();
  topo->setHasChanged(false);

  DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update indexSets start : indexSet0 size : %ld\n", indexSet0->size());
  DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update IndexSets start : indexSet1 size : %ld\n", indexSet1->size());

  // Check indexSet1
  InteractionsGraph::VIterator ui1, ui1end, v1next;
  std::tie(ui1, ui1end) = indexSet1->vertices();

  //Remove interactions from the indexSet1
  for(v1next = ui1; ui1 != ui1end; ui1 = v1next)
  {
    ++v1next;

    SP::Interaction inter1 = indexSet1->bundle(*ui1);
    if(indexSet0->is_vertex(inter1))
    {
      InteractionsGraph::VDescriptor inter1_descr0 = indexSet0->descriptor(inter1);

      assert((indexSet0->color(inter1_descr0) == boost::white_color));

      indexSet0->color(inter1_descr0) = boost::gray_color;
      if(Type::value(*(inter1->nonSmoothLaw())) != Type::EqualityConditionNSL)
      {
        // We assume that the integrator of the ds1 drive the update of the index set
        //SP::OneStepIntegrator Osi = indexSet1->properties(*ui1).osi;
        SP::DynamicalSystem ds1 = indexSet1->properties(*ui1).source;
        OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;

        //if(predictorDeactivate(inter1,i))
        if(osi.removeInteractionFromIndexSet(inter1, i))
        {
          // Interaction is not active
          // ui1 becomes invalid
          indexSet0->color(inter1_descr0) = boost::black_color;

          indexSet1->eraseProperties(*ui1);

          InteractionsGraph::OEIterator oei, oeiend;
          for(std::tie(oei, oeiend) = indexSet1->out_edges(*ui1);
              oei != oeiend; ++oei)
          {
            InteractionsGraph::EDescriptor ed1, ed2;
            std::tie(ed1, ed2) = indexSet1->edges(indexSet1->source(*oei), indexSet1->target(*oei));
            if(ed2 != ed1)
            {
              indexSet1->eraseProperties(ed1);
              indexSet1->eraseProperties(ed2);
            }
            else
            {
              indexSet1->eraseProperties(ed1);
            }
          }
          indexSet1->remove_vertex(inter1);
          /* \warning V.A. 25/05/2012 : Multiplier lambda are only set to zero if they are removed from the IndexSet*/
          inter1->lambda(1)->zero();
          topo->setHasChanged(true);
        }
      }
    }
    else
    {
      // Interaction is not in indexSet0 anymore.
      // ui1 becomes invalid
      indexSet1->eraseProperties(*ui1);
      InteractionsGraph::OEIterator oei, oeiend;
      for(std::tie(oei, oeiend) = indexSet1->out_edges(*ui1);
          oei != oeiend; ++oei)
      {
        InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) = indexSet1->edges(indexSet1->source(*oei), indexSet1->target(*oei));
        if(ed2 != ed1)
        {
          indexSet1->eraseProperties(ed1);
          indexSet1->eraseProperties(ed2);
        }
        else
        {
          indexSet1->eraseProperties(ed1);
        }
      }

      indexSet1->remove_vertex(inter1);
      topo->setHasChanged(true);
    }
  }

  // indexSet0\indexSet1 scan
  InteractionsGraph::VIterator ui0, ui0end;
  //Add interaction in indexSet1
  for(std::tie(ui0, ui0end) = indexSet0->vertices(); ui0 != ui0end; ++ui0)
  {
    if(indexSet0->color(*ui0) == boost::black_color)
    {
      // reset
      indexSet0->color(*ui0) = boost::white_color ;
    }
    else
    {
      if(indexSet0->color(*ui0) == boost::gray_color)
      {
        // reset
        indexSet0->color(*ui0) = boost::white_color;

        assert(indexSet1->is_vertex(indexSet0->bundle(*ui0)));
        /*assert( { !predictorDeactivate(indexSet0->bundle(*ui0),i) ||
          Type::value(*(indexSet0->bundle(*ui0)->nonSmoothLaw())) == Type::EqualityConditionNSL ;
          });*/
      }
      else
      {
        assert(indexSet0->color(*ui0) == boost::white_color);

        SP::Interaction inter0 = indexSet0->bundle(*ui0);
        assert(!indexSet1->is_vertex(inter0));
        bool activate = true;
        if(Type::value(*(inter0->nonSmoothLaw())) != Type::EqualityConditionNSL
            && Type::value(*(inter0->nonSmoothLaw())) != Type::RelayNSL)
        {
          //SP::OneStepIntegrator Osi = indexSet0->properties(*ui0).osi;
          // We assume that the integrator of the ds1 drive the update of the index set
          SP::DynamicalSystem ds1 = indexSet1->properties(*ui0).source;
          OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;

          activate = osi.addInteractionInIndexSet(inter0, i);
        }
        if(activate)
        {
          assert(!indexSet1->is_vertex(inter0));

          // vertex and edges insertion in indexSet1
          indexSet1->copy_vertex(inter0, *indexSet0);
          topo->setHasChanged(true);
          assert(indexSet1->is_vertex(inter0));
        }
      }
    }
  }

  assert(indexSet1->size() <= indexSet0->size());

  DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update indexSets end : indexSet0 size : %ld\n", indexSet0->size());
  DEBUG_PRINTF("TimeStepping::updateIndexSet(unsigned int i). update IndexSets end : indexSet1 size : %ld\n", indexSet1->size());
}

// void TimeStepping::insertNonSmoothProblem(SP::OneStepNSProblem osns)
// {
//   // A the time, a time stepping simulation can only have one non
//   // smooth problem.
//   if((*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY])
//      RuntimeException::selfThrow
//        ("TimeStepping,  insertNonSmoothProblem - A non smooth problem already exist. You can not have more than one.");

//   (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY] = osns;
// }

void TimeStepping::initOSNS()
{
  // === creates links between work vector in OSI and work vector in
  // Interactions
  SP::OneStepIntegrator  osi;

  SP::Topology topo =  _nsds->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);

  InteractionsGraph::VIterator ui, uiend;

  if(!_allNSProblems->empty())  // ie if some Interactions have been
    // declared and a Non smooth problem
    // built.
  {
    //if (_allNSProblems->size()>1)
    //  RuntimeException::selfThrow("TimeStepping::initialize, at the time, a time stepping simulation can not have more than one non smooth problem.");

    // At the time, we consider that for all systems, levelMin is
    // equal to the minimum value of the relative degree - 1 except
    // for degree 0 case where we keep 0.


    // === update all index sets ===
    updateIndexSets();

    // initialization of  OneStepNonSmoothProblem
    for(OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    {
      if(*itOsns)
        (*itOsns)->initialize(shared_from_this());
      else
        RuntimeException::selfThrow("TimeStepping::initOSNS failed. A OneStepNSProblem has not been set. ");
    }
  }
}

void TimeStepping::nextStep()
{
  DEBUG_BEGIN("void TimeStepping::nextStep()\n");
  processEvents();
  DEBUG_END("void TimeStepping::nextStep()\n");
}

void TimeStepping::computeFreeState()
{
  DEBUG_BEGIN("TimeStepping::computeFreeState()\n");
  std::for_each(_allOSI->begin(), _allOSI->end(), std::bind(&OneStepIntegrator::computeFreeState, _1));
  DEBUG_END("TimeStepping::computeFreeState()\n");
}

// compute simulation between current and next event.  Initial
// DS/interaction state is given by memory vectors and final state is
// the one saved in DS/Interaction at the end of this function
void TimeStepping::computeOneStep()
{
  advanceToEvent();
}


void TimeStepping::initializeNewtonLoop()
{
  DEBUG_BEGIN("TimeStepping::initializeNewtonLoop()\n");
  double tkp1 = getTkp1();
  assert(!std::isnan(tkp1));

  for(OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
  {
    (*it)->computeInitialNewtonState();
    (*it)->computeResidu();
  }

  // Predictive contact -- update initial contacts after updating DS positions
  updateWorldFromDS();
  updateInteractions();

  // Changes in updateInteractions may require initialization
  initializeNSDSChangelog();

  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();
  if(indexSet0->size()>0)
  {
    for(OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    {
      (*itOSI)->updateOutput(nextTime());
      (*itOSI)->updateInput(nextTime());
    }
  }

  SP::DynamicalSystemsGraph dsGraph = _nsds->dynamicalSystems();
  for(DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    dsGraph->bundle(*vi)->updatePlugins(tkp1);
  }

  for(OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->computeResidu();

  if(_computeResiduY)
  {
    for(OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    {
      (*itOSI)->computeResiduOutput(tkp1, indexSet0);
    }
  }
  DEBUG_END("TimeStepping::initializeNewtonLoop()\n");
}

void TimeStepping::run()
{
  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of
  // events manager.
  std::cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" <<std::endl;
  while(_eventsManager->hasNextEvent())
  {
    advanceToEvent();

    processEvents();
    count++;
  }
  std::cout << "===== End of " << Type::name(*this) << "simulation. " << count << " events have been processed. ==== " <<std::endl;
}

void TimeStepping::resetLambdas()
{
  if(_resetAllLambda)
  {
    // Initialize lambdas of all interactions.
    SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet(0);
    InteractionsGraph::VIterator ui, uiend, vnext;
    std::tie(ui, uiend) = indexSet0->vertices();
    for(vnext = ui; ui != uiend; ui = vnext)
    {
      ++vnext;
      indexSet0->bundle(*ui)->resetAllLambda();
    }
  }
}

void TimeStepping::advanceToEvent()
{
  DEBUG_PRINTF("TimeStepping::advanceToEvent(). Time =%f\n",getTkp1());
  initialize();
  resetLambdas();
  newtonSolve(_newtonTolerance, _newtonMaxIteration);
}

/*update of the nabla */
/*discretisation of the Interactions */
void   TimeStepping::prepareNewtonIteration()
{
  DEBUG_BEGIN("TimeStepping::prepareNewtonIteration()\n");
  double tkp1 = getTkp1();
  for(OSIIterator itosi = _allOSI->begin();
      itosi != _allOSI->end(); ++itosi)
  {
    (*itosi)->prepareNewtonIteration(tkp1);
  }

  DEBUG_END("TimeStepping::prepareNewtonIteration()\n");
}

void TimeStepping::saveYandLambdaInOldVariables()
{
  // Temp FP : saveInOldVar was called for each osns and each osns call
  // swapInOldVar for all interactions in the nsds.
  // ==> let's do it only once, by the simu.

  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();
  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    //indexSet0->bundle(*ui)->swapInMemory();
    indexSet0->bundle(*ui)->swapInOldVariables();;
  }
}

void TimeStepping::displayNewtonConvergenceInTheLoop()
{
  if(_displayNewtonConvergence)
  {
    std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonNbIterations =" << _newtonNbIterations << std::endl;
    std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonResiduDSMax =" << _newtonResiduDSMax << std::endl;
    if(_computeResiduY)
    {
      std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonResiduYMax =" << _newtonResiduYMax << std::endl;
    }
    else
    {
      std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonResiduYMax =" << "not computed" << std::endl;
    }
    if(_computeResiduR)
    {
      std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonResiduRMax =" << _newtonResiduRMax << std::endl;
    }
    else
    {
      std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonResiduRMax =" << "not computed" << std::endl;
    }
  }
  else
  {
    DEBUG_PRINTF("# _newtonNbIterations = %i\n",_newtonNbIterations);
    DEBUG_PRINTF("# _newtonResiduDSMax = %12.8e\t",_newtonResiduDSMax);
    DEBUG_PRINTF("# _newtonResiduYMax = %12.8e\t",_newtonResiduYMax);
    DEBUG_PRINTF("# _newtonResiduRMax = %12.8e\n",_newtonResiduRMax);
  }
}
void TimeStepping::displayNewtonConvergenceAtTheEnd(int info, unsigned int maxStep)
{
  if(_displayNewtonConvergence)
  {
    std::cout << "[kernel] TimeStepping::newtonSolve --  _newtonCumulativeNbIterations =" << _newtonCumulativeNbIterations << std::endl;
  }
  else
  {
    DEBUG_PRINTF("# _newtonCumulativeNbIterations= %i\n",_newtonCumulativeNbIterations);
  }

  if(!_isNewtonConverge)
  {
    if(_warnOnNonConvergence)
      std::cout << "[kernel] TimeStepping::newtonSolve -- Newton process stopped:"
                << "max. number of steps ("
                << maxStep
                <<") reached at accuracy = "
                << _newtonResiduDSMax
                <<"."
                << std::endl ;
    if(info && _warnOnNonConvergence)
      std::cout << "[kernel] TimeStepping::newtonSolve -- nonsmooth solver failed." <<std::endl ;
  }
}
void TimeStepping::newtonSolve(double criterion, unsigned int maxStep)
{

  DEBUG_BEGIN("TimeStepping::newtonSolve(double criterion, unsigned int maxStep)\n");
  _isNewtonConverge = false;
  _newtonNbIterations = 0; // number of Newton iterations
  int info = 0;
  bool isLinear  = _nsds->isLinear();

  initializeNewtonLoop();

  if((_newtonOptions == SICONOS_TS_LINEAR || _newtonOptions == SICONOS_TS_LINEAR_IMPLICIT)
      || isLinear)
  {
    _newtonNbIterations++;
    DEBUG_PRINTF("TimeStepping::newtonSolve(). _newtonNbIterations = %i\n", _newtonNbIterations);
    prepareNewtonIteration();
    computeFreeState();
    //bool hasNSProblems = (!_allNSProblems->empty() &&   indexSet0.size() > 0) ? true : false;
    bool hasNSProblems = (!_allNSProblems->empty()) ? true : false;
    if(hasNSProblems)
      info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
    // Check output from solver (convergence or not ...)
    if(!checkSolverOutput)
      DefaultCheckSolverOutput(info);
    else
      checkSolverOutput(info, this);

    update();

    hasNSProblems = (!_allNSProblems->empty()) ? true : false;
    if(hasNSProblems)
      saveYandLambdaInOldVariables();
  }

  else if(_newtonOptions == SICONOS_TS_NONLINEAR)
  {
    //  while((!_isNewtonConverge)&&(_newtonNbIterations < maxStep)&&(!info))
    //_isNewtonConverge = newtonCheckConvergence(criterion);
    while((!_isNewtonConverge) && (_newtonNbIterations < maxStep))
    {
      _newtonNbIterations++;
      info=0;

      prepareNewtonIteration();
      computeFreeState();

      if(info && _warnOnNonConvergence)
        std::cout << "New Newton loop because of nonsmooth solver failed\n" <<std::endl;

      // if there is not any Interaction at
      // the beginning of the simulation _allNSProblems may not be
      // empty here (check with SpaceFilter and one disk not on
      // the ground : MultiBodyTest::t2)

      // if((*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->simulation())
      // is also relevant here.
      //InteractionsGraph& indexSet0 = *_nsds->topology()->indexSet0();
      //bool hasNSProblems = (!_allNSProblems->empty() &&   indexSet0.size() > 0) ? true : false;

      bool hasNSProblems = (!_allNSProblems->empty()) ? true : false;
      if(hasNSProblems)
      {
        info = computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);
        // Check output from solver (convergence or not ...)
        if(!checkSolverOutput)
          DefaultCheckSolverOutput(info);
        else
          checkSolverOutput(info, this);
      }

      updateInput();
      updateState();

      if(!_isNewtonConverge && _newtonNbIterations < maxStep)
      {
        updateOutput();
      }

      _isNewtonConverge = newtonCheckConvergence(criterion);

      if(!_isNewtonConverge && !info)
      {
        if(hasNSProblems)
          saveYandLambdaInOldVariables();
      }
      displayNewtonConvergenceInTheLoop();
    } // End of the Newton Loop

    _newtonCumulativeNbIterations += _newtonNbIterations;

    displayNewtonConvergenceAtTheEnd(info, maxStep);
  }
  else
    RuntimeException::selfThrow("TimeStepping::NewtonSolve failed. Unknown newtonOptions: " + std::to_string(_newtonOptions));
  DEBUG_END("TimeStepping::newtonSolve(double criterion, unsigned int maxStep)\n");
}

bool TimeStepping::newtonCheckConvergence(double criterion)
{
  bool checkConvergence = true;
  //_relativeConvergenceCriterionHeld is true means that the RCC is
  //activated, and the relative criteron helds.  In this case the
  //newtonCheckConvergence has to return true. End of the Newton
  //iterations
  if(_relativeConvergenceCriterionHeld)
  {
    return true;
  }
  // get the nsds indicator of convergence
  // We compute cvg = abs(xi+1 - xi)/xi and if cvg < criterion
  //  if (nsdsConverge < criterion )

  double residu = 0.0;
  _newtonResiduDSMax = 0.0;
  for(OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
  {
    residu = (*it)->computeResidu();

    if(residu > _newtonResiduDSMax) _newtonResiduDSMax = residu;
    if(residu > criterion)
    {
      checkConvergence = false;
      //break;
    }
  }




  if(_computeResiduY)
  {
    //check residuy.
    _newtonResiduYMax = 0.0;
    residu = 0.0;

    SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();
    for(OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    {
      residu = std::max(residu,(*itOSI)->computeResiduOutput(getTkp1(), indexSet0));
    }


//     InteractionsGraph::VIterator ui, uiend;
//     SP::Interaction inter;
//     for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
//     {
//       inter = indexSet0->bundle(*ui);
//       VectorOfVectors& workV = *indexSet0->properties(*ui).workVectors;

//       inter->computeResiduY(, workV);
//     inter->residuY()->norm2();
    if(residu > _newtonResiduYMax) _newtonResiduYMax = residu;
    if(residu > criterion)
      checkConvergence = false;
  }

  if(_computeResiduR)
  {
    //check residur.
    _newtonResiduRMax = 0.0;
    residu = 0.0;
    SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();

    for(OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    {
      residu = std::max(residu,(*itOSI)->computeResiduInput(getTkp1(), indexSet0));
    }


    // InteractionsGraph::VIterator ui, uiend;
    // SP::Interaction inter;
    // for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    // {
    //   inter = indexSet0->bundle(*ui);
    //   VectorOfBlockVectors& DSlink = *indexSet0->properties(*ui).DSlink;
    //   VectorOfVectors& workV = *indexSet0->properties(*ui).workVectors;

    //   inter->computeResiduR(getTkp1(), DSlink, workV);
    //   // TODO support other DS
    if(residu > _newtonResiduRMax) _newtonResiduRMax = residu;
    if(residu > criterion)
    {
      checkConvergence = false;
    }
  }


  return(checkConvergence);
}

void TimeStepping::DefaultCheckSolverOutput(int info)
{
  // info = 0 => ok
  // else: depend on solver
  if(info != 0)
  {
    std::cout << "TimeStepping::DefaultCheckSolverOutput:" << std::endl;
    std::cout << "Non smooth solver warning : output message from numerics solver is equal to " << info << std::endl;
    //       std::cout << "=> may have failed? (See Numerics solver documentation for details on the message meaning)." <<std::endl;
    //      std::cout << "=> may have failed? (See Numerics solver documentation for details on the message meaning)." <<std::endl;
    //     RuntimeException::selfThrow(" Non smooth problem, solver convergence failed ");
    /*      if(info == 1)
            std::cout <<" reach max iterations number with solver " << solverName <<std::endl;
            else if (info == 2)
            {
            if (solverName == "LexicoLemke" || solverName == "CPG" || solverName == "NLGS")
            RuntimeException::selfThrow(" negative diagonal term with solver "+solverName);
            else if (solverName == "QP" || solverName == "NSQP" )
            RuntimeException::selfThrow(" can not satisfy convergence criteria for solver "+solverName);
            else if (solverName == "Latin")
            RuntimeException::selfThrow(" Choleski factorisation failed with solver Latin");
            }
            else if (info == 3 && solverName == "CPG")
            std::cout << "pWp null in solver CPG" <<std::endl;
            else if (info == 3 && solverName == "Latin")
            RuntimeException::selfThrow("Null diagonal term with solver Latin");
            else if (info == 5 && (solverName == "QP" || solverName == "NSQP"))
            RuntimeException::selfThrow("Length of working array insufficient in solver "+solverName);
            else
            RuntimeException::selfThrow("Unknown error type in solver "+ solverName);
    */
  }
}

void TimeStepping::setCheckSolverFunction(CheckSolverFPtr newF)
{
  checkSolverOutput = newF;
}
