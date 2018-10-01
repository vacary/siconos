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
#include "LagrangianDS.hpp"
#include "BlockVector.hpp"
#include "BlockMatrix.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
#include <iostream>

void LagrangianDS::_init(SP::SiconosVector position, SP::SiconosVector velocity)
{
  assert(_ndof > 0 && "lagrangian dynamical system dimension should be greater than 0.");

  // Set initial conditions
  _q0 = position;
  _velocity0 = velocity;

  // -- Memory allocation for vector and matrix members --
  _q.resize(3);
  _q[0].reset(new SiconosVector(*_q0));
  _q[1].reset(new SiconosVector(*_velocity0));

  /** \todo lazy Memory allocation */
  _p.resize(3);
  _p[1].reset(new SiconosVector(_ndof));

  _zeroPlugin();
}


// Build from initial state only
LagrangianDS::LagrangianDS(SP::SiconosVector q0, SP::SiconosVector v0):
  DynamicalSystem(2 * q0->size()), _ndof(v0->size()),
  _hasConstantMass(true), _hasConstantK(false), _hasConstantC(false), _hasConstantFExt(true)
{
  // Initial conditions
  _init(q0, v0);
}

// From initial state and constant mass matrix, \f$ M\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosMatrix newMass):
  DynamicalSystem(2 * q0->size()), _ndof(v0->size()),
  _hasConstantMass(true),  _hasConstantK(false), _hasConstantC(false), _hasConstantFExt(true)

{
  _init(q0, v0);
  // Mass matrix
  _mass = newMass;
}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector q0, SP::SiconosVector v0, const std::string& massName):
  DynamicalSystem(), _ndof(q0->size()),
  _hasConstantMass(false), _hasConstantK(false), _hasConstantC(false), _hasConstantFExt(true)
{
  _init(q0, v0);
  // Mass
  _mass.reset(new SimpleMatrix(_ndof, _ndof));
  setComputeMassFunction(SSLH::getPluginName(massName), SSLH::getPluginFunctionName(massName));
}

void LagrangianDS::_zeroPlugin()
{
  _pluginMass.reset(new PluggedObject());
  _pluginFInt.reset(new PluggedObject());
  _pluginFExt.reset(new PluggedObject());
  _pluginFGyr.reset(new PluggedObject());
  _pluginJacqFInt.reset(new PluggedObject());
  _pluginJacqDotFInt.reset(new PluggedObject());
  _pluginJacqFGyr.reset(new PluggedObject());
  _pluginJacqDotFGyr.reset(new PluggedObject());
}

void LagrangianDS::initializeNonSmoothInput(unsigned int level)
{
  if(!_p[level])
    _p[level].reset(new SiconosVector(_ndof));
}

void LagrangianDS::resetToInitialState()
{
  if(_q0)
  {
    *(_q[0]) = *_q0;
  }
  else
    RuntimeException::selfThrow("LagrangianDS::resetToInitialState - initial position _q0 is null");
  if(_velocity0)
  {
    *(_q[1]) = *_velocity0;
  }
  else
    RuntimeException::selfThrow("LagrangianDS::resetToInitialState - initial velocity _velocity0 is null");
}

void LagrangianDS::init_generalized_coordinates(unsigned int level)
{
  assert(level>1);
  if(!_q[level])
    _q[level].reset(new SiconosVector(_ndof));
}


void LagrangianDS::init_inverse_mass()
{
  if(_mass && !_inverseMass)
  {
    computeMass();
    _inverseMass.reset(new SimpleMatrix(*_mass));
  }
}

void LagrangianDS::update_inverse_mass()
{
  if(_mass && _inverseMass && !_hasConstantMass)
  {
    computeMass();
    *_inverseMass = *_mass;
  }
}

void LagrangianDS::init_forces()
{
  DEBUG_BEGIN("LagrangianDS::init_forces()\n");
  // Allocate memory for forces and its jacobians.
  if(_fInt || _fExt || _fGyr)
  {
    if(!_forces)
      _forces.reset(new SiconosVector(_ndof));
  }

  if(_fInt || _fGyr)
  {
    if(!_jacobianqForces)
    {
      /* Warning: we are forcing the dense storage */
      _jacobianqForces.reset(new SimpleMatrix(_ndof, _ndof));
      DEBUG_PRINT("_jacobianqForces.reset(new SimpleMatrix(_ndof, _ndof));\n");
      DEBUG_EXPR(_jacobianqForces->display(););
    }
    if(!_jacobianqDotForces)
      _jacobianqDotForces.reset(new SimpleMatrix(_ndof, _ndof));
  }
  DEBUG_END("LagrangianDS::init_forces()\n");
}

void LagrangianDS::initRhs(double time)
{
  DEBUG_BEGIN("LagrangianDS::initRhs(double time)\n");
  // dim
  _n = 2 * _ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianfx.

  // Initial conditions and state

  // WARNING : this function is supposed to be called
  // by the OneStepIntegrator, and maybe several times for the same DS
  // if the system is involved in more than one interaction. So, we must check
  // if p2 and q2 already exist to be sure that DSlink won't be lost.

  _x0.reset(new SiconosVector(*_q0, *_velocity0));

  _x[0].reset(new SiconosVector(*_q[0], *_q[1]));

  if(!_q[2])
    _q[2].reset(new SiconosVector(_ndof));

  _x[1].reset(new SiconosVector(*_q[1], *_q[2]));

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
  _rhsMatrices.resize(numberOfRhsMatrices);

  if(!_p[2])
    _p[2].reset(new SiconosVector(_ndof));

  init_forces();
  init_inverse_mass();

  computeRhs(time);

  bool flag1 = false, flag2 = false;
  if(_jacobianqForces)
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqForces(time);

    _rhsMatrices[jacobianXBloc10].reset(new SimpleMatrix(*_jacobianqForces));
    _inverseMass->PLUForwardBackwardInPlace(*_rhsMatrices[jacobianXBloc10]);
    flag1 = true;
  }

  if(_jacobianqDotForces)
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianqDotForces(time);
    _rhsMatrices[jacobianXBloc11].reset(new SimpleMatrix(*_jacobianqDotForces));
    _inverseMass->PLUForwardBackwardInPlace(*_rhsMatrices[jacobianXBloc11]);
    flag2 = true;
  }

  if(!_rhsMatrices[zeroMatrix])
    _rhsMatrices[zeroMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::ZERO));
  if(!_rhsMatrices[idMatrix])
    _rhsMatrices[idMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::IDENTITY));

  if(flag1 && flag2)
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[jacobianXBloc10], _rhsMatrices[jacobianXBloc11]));
  else if(flag1)  // flag2 = false
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[jacobianXBloc10], _rhsMatrices[zeroMatrix]));
  else if(flag2)  // flag1 = false
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[zeroMatrix], _rhsMatrices[jacobianXBloc11]));
  else
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[zeroMatrix], _rhsMatrices[zeroMatrix]));
  DEBUG_EXPR(display(););
  DEBUG_END("LagrangianDS::initRhs(double time)\n");
}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if(! _q[0])
    _q[0].reset(new SiconosVector(newValue));
  else
    *_q[0] = newValue;
}

void LagrangianDS::setQPtr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");
  _q[0] = newPtr;

}

void LagrangianDS::setQ0(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if(! _q0)
    _q0.reset(new SiconosVector(newValue));
  else
    *_q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");
  _q0 = newPtr;
}

void LagrangianDS::setVelocity0(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if(! _velocity0)
    _velocity0.reset(new SiconosVector(newValue));
  else
    *_velocity0 = newValue;
}

void LagrangianDS::setVelocity(const SiconosVector& newValue)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if(! _q[1])
    _q[1].reset(new SiconosVector(newValue));
  else
    *_q[1] = newValue;
}

void LagrangianDS::setVelocityPtr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");
  _q[1] = newPtr;
}


void LagrangianDS::setVelocity0Ptr(SP::SiconosVector newPtr)
{
  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");
  _velocity0 = newPtr;
}

void LagrangianDS::setP(const SiconosVector& newValue, unsigned int level)
{
  if(newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if(! _p[level])
    _p[level].reset(new SiconosVector(newValue));
  else
    *(_p[level]) = newValue;
}

void LagrangianDS::setPPtr(SP::SiconosVector newPtr, unsigned int level)
{

  if(newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");
  _p[level] = newPtr;
}

void LagrangianDS::computeMass()
{
  DEBUG_BEGIN("LagrangianDS::computeMass()\n");
  DEBUG_EXPR(_q[0]->display());
  if(!_hasConstantMass)
  {
    if(_mass && _pluginMass->fPtr)
    {
      ((FPtr7)_pluginMass->fPtr)(_ndof, &(*_q[0])(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
      _mass->resetLU();
    }
  }
  DEBUG_EXPR(_mass->display());
  DEBUG_END("LagrangianDS::computeMass()\n");
}

void LagrangianDS::computeMass(SP::SiconosVector position)
{
  if(_mass && !_hasConstantMass && _pluginMass->fPtr)
  {
    ((FPtr7)_pluginMass->fPtr)(_ndof, &(*position)(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
    _mass->resetLU();
  }
}

/** This function has been added to avoid Swig director to wrap _fExt into
 * numpy.array when we call
 * LagrangianDS::computeFExt(double time, SP::SiconosVector fExt)
 * that calls in turn computeFExt(time, _fExt);
 */
static
void computeFExt_internal(double time, bool hasConstantFExt,
                          unsigned int ndof,
                          SP::PluggedObject pluginFExt, SP::SiconosVector fExt_attributes,
                          SP::SiconosVector fExt,
                          SP::SiconosVector z)
{
/* if the pointer has been set to an external vector
 * after setting the plugin, we do not call the plugin */
  if(hasConstantFExt)
  {
    if(fExt != fExt_attributes)
      *fExt = *fExt_attributes;
  }
  else if(pluginFExt->fPtr)
    ((VectorFunctionOfTime)pluginFExt->fPtr)(time, ndof, &(*fExt)(0), z->size(), &(*z)(0));

}

void LagrangianDS::computeFExt(double time)
{
  computeFExt_internal(time,_hasConstantFExt,
                       _ndof,
                       _pluginFExt, _fExt, _fExt, _z);
}

void LagrangianDS::computeFExt(double time, SP::SiconosVector fExt)
{
  computeFExt_internal(time,_hasConstantFExt,
                       _ndof,
                       _pluginFExt, _fExt, fExt, _z);
}

/** This function has been added to avoid Swig director to wrap _FInt into
* numpy.array  when we call
*  LagrangianDS::computeFInt(double time, SP::SiconosVector fInt)
* that calls in turn computeFInt(time, ...;
*/
static
void computeFInt_internal(double time, bool hasConstantK, bool hasConstantC,
                          unsigned int ndof,
                          SP::SiconosVector q, SP::SiconosVector v,
                          SP::SiconosMatrix K, SP::SiconosMatrix C,
                          SP::PluggedObject pluginFInt, SP::SiconosVector fInt_attributes,
                          SP::SiconosVector fInt,
                          SP::SiconosVector z)
{
  /* if the pointer has been set to an external matrix for K or C
   * after setting the plugin, we do not call the plugin */
  if(hasConstantK || hasConstantC)
  {
    if(fInt != fInt_attributes)
      *fInt = *fInt_attributes;
    fInt->zero();
    if (hasConstantK)
      prod(1.0, *K , *q, *fInt, false);
    if (hasConstantC)
      prod(1.0, *C , *v, *fInt, false);
  }
  else if(pluginFInt->fPtr)
    ((FPtr6)pluginFInt->fPtr)(time, ndof, &(*q)(0), &(*v)(0), &(*fInt)(0), z->size(), &(*z)(0));
}

void LagrangianDS::computeFInt(double time)
{
  computeFInt_internal(time, _hasConstantK, _hasConstantC,
                       _ndof,
                       _q[0], _q[1],
                       _K, _C,
                       _pluginFInt, _fInt, _fInt, _z);
}

void LagrangianDS::computeFInt(double time,
                               SP::SiconosVector position,
                               SP::SiconosVector velocity,
                               SP::SiconosVector fInt)
{
  computeFInt_internal(time, _hasConstantK, _hasConstantC,
                       _ndof,
                       position, velocity,
                       _K,_C,
                       _pluginFInt, _fInt, fInt, _z);
}




void LagrangianDS::computeFGyr()
{
  if(_fGyr && _pluginFGyr->fPtr)
    ((FPtr5)_pluginFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_fGyr)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFGyr(SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_fGyr && _pluginFGyr->fPtr)
    ((FPtr5)_pluginFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0), &(*_fGyr)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::setK(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != _ndof || newValue.size(1) != _ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (!_K)
    _K.reset(new SimpleMatrix(newValue));
  else
    *_K = newValue;
}

void LagrangianDS::setC(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != _ndof || newValue.size(1) != _ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setC: inconsistent input matrix size ");

  if (!_C)
    _C.reset(new SimpleMatrix(newValue));
  else
    *_C = newValue;
}

void LagrangianDS::computeJacobianFIntq(double time)
{
  DEBUG_BEGIN("LagrangianDS::computeJacobianFIntq()\n");
  DEBUG_EXPR(_q[0]->display());
  DEBUG_EXPR(_q[1]->display());
  if (!_hasConstantK)
  {
    if(_K&& _pluginJacqFInt->fPtr)
      ((FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_K)(0, 0), _z->size(), &(*_z)(0));
  }
  DEBUG_EXPR(if(_K) _K->display(););
  DEBUG_END("LagrangianDS::computeJacobianFIntq()\n");
}

void LagrangianDS::computeJacobianFIntqDot(double time)
{
  DEBUG_BEGIN("LagrangianDS::computeJacobianFIntqDot()\n");
  DEBUG_EXPR(_q[0]->display());
  DEBUG_EXPR(_q[1]->display());
  DEBUG_EXPR(_z->display());
  if(!_hasConstantC)
  {
    if(_C && _pluginJacqDotFInt->fPtr)
      ((FPtr6)_pluginJacqDotFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_C)(0, 0), _z->size(), &(*_z)(0));
  }
  DEBUG_EXPR(if(_C) _C->display(););
  DEBUG_END("LagrangianDS::computeJacobianFIntqDot()\n");
}


void LagrangianDS::computeJacobianFIntq(double time, SP::SiconosVector position, SP::SiconosVector velocity)
{
  DEBUG_BEGIN("LagrangianDS::computeJacobianFIntq()\n");
  DEBUG_EXPR(position->display());
  DEBUG_EXPR(velocity->display());
  if (!_hasConstantK)
  {
    if(_K && _pluginJacqFInt->fPtr)
      ((FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*position)(0), &(*velocity)(0), &(*_K)(0, 0), _z->size(), &(*_z)(0));
  }
  DEBUG_EXPR(if(_K) _K->display(););
  DEBUG_END("LagrangianDS::computeJacobianFIntq()\n");
}
void LagrangianDS::computeJacobianFIntqDot(double time, SP::SiconosVector position, SP::SiconosVector velocity)
{
  if (!_hasConstantC)
  {
    if(_C && _pluginJacqDotFInt->fPtr)
      ((FPtr6)_pluginJacqDotFInt->fPtr)(time, _ndof, &(*position)(0), &(*velocity)(0), &(*_C)(0, 0), _z->size(), &(*_z)(0));
  }
}

void LagrangianDS::computeJacobianFGyrq()
{
  if(_pluginJacqFGyr->fPtr)
    ((FPtr5)_pluginJacqFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFGyrq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianFGyrqDot()
{
  if(_jacobianFGyrqDot && _pluginJacqDotFGyr->fPtr)
    ((FPtr5)_pluginJacqDotFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFGyrqDot)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFGyrq(SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_jacobianFGyrq && _pluginJacqFGyr->fPtr)
    ((FPtr5)_pluginJacqFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFGyrq)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFGyrqDot(SP::SiconosVector position, SP::SiconosVector velocity)
{
  if(_jacobianFGyrqDot && _pluginJacqDotFGyr->fPtr)
    ((FPtr5)_pluginJacqDotFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFGyrqDot)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeRhs(double time)
{
  DEBUG_BEGIN("LagrangianDS::computeRhs(double time)");
  *_q[2] = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  // if(_forces)
  //   {
  computeForces(time, _q[0], _q[1]);
  *_q[2] += *_forces;
  DEBUG_EXPR(_forces->display(););
  //#  }

  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during initialization and saved into _inverseMass
  // -- Case 2: mass is not constant, we copy it into _inverseMass
  // Then we proceed with PLUForwardBackward.
  // mass and inv(mass) computation
  if(_mass && !_hasConstantMass)  // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
  {
    computeMass();
    *_inverseMass = *_mass;
  }

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianDS children
  if(_inverseMass)
    _inverseMass->PLUForwardBackwardInPlace(*_q[2]);

  _x[1]->setBlock(0, *_q[1]);
  _x[1]->setBlock(_ndof, *_q[2]);
  DEBUG_END("LagrangianDS::computeRhs(double time)");
}

void LagrangianDS::computeJacobianRhsx(double time)
{
  if(!_hasConstantMass)
    computeMass();

  //  if(mass->isPlugged()) : mass may b not plugged in LagrangianDS children

  if(_jacobianqForces || _jacobianqDotForces)
  {
    if(!_hasConstantMass) // else inverseMass is already uptodate
      *_inverseMass = *_mass;
  }

  if(_jacobianqForces)
  {
    /** \warning the Jacobian of the inverse of the mass matrix
     * w.r.t q is not taken into account */

    SP::SiconosMatrix bloc10 = _jacxRhs->block(1, 0);
    computeJacobianqForces(time);
    *bloc10 = *_jacobianqForces;
    _inverseMass->PLUForwardBackwardInPlace(*bloc10);
  }

  if(_jacobianqDotForces)
  {
    SP::SiconosMatrix bloc11 = _jacxRhs->block(1, 1);
    computeJacobianqDotForces(time);
    *bloc11 = *_jacobianqDotForces;
    _inverseMass->PLUForwardBackwardInPlace(*bloc11);
  }
}

void LagrangianDS::computeForces(double time,
                                 SP::SiconosVector position,
                                 SP::SiconosVector velocity)
{
  if(!_forces)
  {
    _forces.reset(new SiconosVector(_ndof));
  }
  else
    _forces->zero();

  // 1 - Computes the required function
  computeFInt(time, position, velocity, _fInt);
  computeFExt(time);
  computeFGyr(position, velocity);

  // seems ok.
  // if (_forces.use_count() == 1)
  // {
  //   //if not that means that fL is already (pointer-)connected with
  //   // either fInt, FGyr OR fExt.
  //_forces->zero();

  if(_fInt)
    *_forces -= *_fInt;

  if(_fExt)
    *_forces += *_fExt;

  if(_fGyr)
    *_forces -= *_fGyr;
  // }

}

void LagrangianDS::computeJacobianqForces(double time)
{
  DEBUG_BEGIN("LagrangianDS::computeJacobianqForces(double time)\n");
  if (!_hasConstantK)
  {
    if(_jacobianqForces)
    {
      computeJacobianFIntq(time);
      computeJacobianFGyrq();

      // not true!
      // if( jacobianFL[i].use_count() == 1 )
      {
        //if not that means that jacobianFL[i] is already (pointer-)connected with
        // either jacobianFInt or jacobianFGyr
        _jacobianqForces->zero();
        if(_K)
          *_jacobianqForces -= *_K;
        if(_jacobianFGyrq)
          *_jacobianqForces -= *_jacobianFGyrq;
      }
    }
  }
  //else nothing.
  DEBUG_END("LagrangianDS::computeJacobianqForces(double time)\n");
}
void LagrangianDS::computeJacobianqDotForces(double time)
{
  if(_jacobianqDotForces)
  {
    computeJacobianFIntqDot(time);
    computeJacobianFGyrqDot();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianFGyr
      _jacobianqDotForces->zero();
      if(_C)
        *_jacobianqDotForces -= *_C;
      if(_jacobianFGyrqDot)
        *_jacobianqDotForces -= *_jacobianFGyrqDot;
    }
  }
  //else nothing.
}
// void LagrangianDS::computeJacobianZFL( double time){
//    RuntimeException::selfThrow("LagrangianDS::computeJacobianZFL - not implemented");
// }


// --- Functions for memory handling ---
void LagrangianDS::initMemory(unsigned int steps)
{
  DEBUG_PRINTF("LagrangianDS::initMemory(unsigned int steps) with steps = %i\n", steps);
  if(steps == 0)
    std::cout << "Warning : LagragianDS::initMemory with size equal to zero" <<std::endl;
  else
  {
    _qMemory.setMemorySize(steps, _ndof);
    _velocityMemory.setMemorySize(steps, _ndof);
    _forcesMemory.setMemorySize(steps, _ndof);
    _pMemory.resize(3);

    //TODO : initMemory in graph + f(OSI/level)
    for(unsigned int level=0; level<3; ++level)
    {
      if(_pMemory[level].size()==0)
        _pMemory[level].setMemorySize(steps, _ndof);
    }

    //swapInMemory();
  }
}

void LagrangianDS::swapInMemory()
{
  _qMemory.swap(*_q[0]);
  _velocityMemory.swap(*_q[1]);
  if (_forces)
    _forcesMemory.swap(*_forces);

  // initialization of the reaction force due to the non smooth law
  // note: these are a no-op if either memory or vector is null
  _pMemory[0].swap(_p[0]);
  _pMemory[1].swap(_p[1]);
  _pMemory[2].swap(_p[2]);
  _xMemory.swap(_x[0]);
}

void LagrangianDS::resetAllNonSmoothParts()
{
  if(_p[0])
    _p[0]->zero();
  if(_p[1])
    _p[1]->zero();
  if(_p[2])
    _p[2]->zero();
}

void LagrangianDS::resetNonSmoothPart(unsigned int level)
{
  if(level < LEVELMAX)
    if(_p[level])
      _p[level]->zero();
}

void LagrangianDS::computePostImpactVelocity()
{
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  DEBUG_BEGIN("LagrangianDS::computePostImpactV()\n");
  SiconosVector tmp(*_p[1]);
  if(_inverseMass)
    _inverseMass->PLUForwardBackwardInPlace(tmp);
  *_q[1] += tmp;  // v+ = v- + p
  DEBUG_BEGIN("LagrangianDS::computePostImpactV() END \n");
}

void LagrangianDS::setComputeFGyrFunction(const std::string& pluginPath, const std::string&  functionName)
{
  _pluginFGyr->setComputeFunction(pluginPath, functionName);
  if(!_fGyr)
    _fGyr.reset(new SiconosVector(_ndof));
  init_forces();
}

void LagrangianDS::setComputeFGyrFunction(FPtr5 fct)
{
  _pluginFGyr->setComputeFunction((void *)fct);
  if(!_fGyr)
    _fGyr.reset(new SiconosVector(_ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
  if(!_K)
    _K.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqDotFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqDotFInt->setComputeFunction(pluginPath, functionName);
  if(!_C)
    _C.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqFunction(FPtr6 fct)
{
  _pluginJacqFInt->setComputeFunction((void *)fct);
  if(!_K)
    _K.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFIntqDotFunction(FPtr6 fct)
{
  _pluginJacqDotFInt->setComputeFunction((void *)fct);
  if(!_C)
    _C.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqFGyr->setComputeFunction(pluginPath, functionName);
  if(!_jacobianFGyrq)
    _jacobianFGyrq.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqDotFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqDotFGyr->setComputeFunction(pluginPath, functionName);
  if(!_jacobianFGyrqDot)
    _jacobianFGyrqDot.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqFunction(FPtr5 fct)
{
  _pluginJacqFGyr->setComputeFunction((void *)fct);
  if(!_jacobianFGyrq)
    _jacobianFGyrq.reset(new SimpleMatrix(_ndof, _ndof));
  init_forces();
}

void LagrangianDS::setComputeJacobianFGyrqDotFunction(FPtr5 fct)
{
  _pluginJacqDotFGyr->setComputeFunction((void *)fct);
  if(!_jacobianFGyrqDot)
    _jacobianFGyrqDot.reset(new SimpleMatrix(_ndof, _ndof));
}

double LagrangianDS::computeKineticEnergy()
{
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  SP::SiconosVector velo = velocity();
  assert(velo);
  DEBUG_EXPR(velo->display());

  SP::SiconosVector tmp(new SiconosVector(*velo));
  if(_mass)
    prod(*_mass, *velo, *tmp, true);

  double K =0.5*inner_prod(*tmp,*velo);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("LagrangianDS::computeKineticEnergy()\n");
  return K;
}

void LagrangianDS::setBoundaryConditions(SP::BoundaryCondition newbd)
{
  if(!_boundaryConditions)
  {
    std::cout << "Warning : LagrangianDS::setBoundaryConditions. old boundary conditions were pre-existing" <<std::endl;
  }
  _boundaryConditions = newbd;
  _reactionToBoundaryConditions.reset(new SiconosVector(_boundaryConditions->velocityIndices()->size()));
};


void LagrangianDS::display(bool brief) const
{
  std::cout << "=====> Lagrangian System display (number: " << _number << ")." <<std::endl;
  std::cout << "- _ndof : " << _ndof <<std::endl;
  std::cout << "- q " <<std::endl;
  if(_q[0]) _q[0]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- q0 " <<std::endl;
  if(_q0) _q0->display();
  std::cout << "- velocity " <<std::endl;
  if(_q[1]) _q[1]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- acceleration " <<std::endl;
  if(_q[2]) _q[2]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- v0 " <<std::endl;
  if(_velocity0) _velocity0->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[0] " <<std::endl;
  if(_p[0]) _p[0]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[1] " <<std::endl;
  if(_p[1]) _p[1]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[2] " <<std::endl;
  if(_p[2]) _p[2]->display();
  else std::cout << "-> NULL" <<std::endl;

  if (!brief)
  {
    std::cout << "- fExt " <<std::endl;
    if(_fExt) _fExt->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "- fInt " <<std::endl;
    if(_fInt) _fInt->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "- K " <<std::endl;
    if(_K) _K->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "- C " <<std::endl;
    if(_C) _C->display();
    else std::cout << "-> NULL" <<std::endl;

  }
  std::cout << "===================================== " <<std::endl;
}
