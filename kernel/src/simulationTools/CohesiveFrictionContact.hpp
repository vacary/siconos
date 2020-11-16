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
/*! \file
  Fricton-Contact Non-Smooth Problem
*/
#ifndef CohesiveFrictionContact_H
#define CohesiveFrictionContact_H

#include "FrictionContact.hpp"

#include <FrictionContactProblem.h>
#include <Friction_cst.h>
#include <fc2d_Solvers.h>
#include <fc3d_Solvers.h>

/** Pointer to function of the type used for drivers for CohesiveFrictionContact problems in Numerics */
// typedef int (*Driver)(FrictionContactProblem*, double*, double*, SolverOptions*);
// TYPEDEF_SPTR(FrictionContactProblem)


/** Formalization and Resolution of a Friction-Contact Problem

  This class is devoted to the formalization and the resolution of
  friction contact problems defined by :


  \rst

  .. math::

     velocity =  q + M reaction \\
     \\
     velocity \geq 0, reaction \geq 0,  reaction^{T} velocity =0

  \endrst

  and a Coulomb friction law.

  With:
     - \f$velocity \in R^{n} \f$  and \f$reaction \in R^{n} \f$ the unknowns,
     - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$

  The dimension of the problem (2D or 3D) is given by the variable contactProblemDim and the proper
  Numerics driver will be called according to this value.

  \b Construction: just set Numerics Solver id

  Main functions:

  \b Usage:
  - compute(time) formalize, solve and post-process the problem.

  pre- and post-pro are common to all LinearOSNS and defined in this class.


 */
class CohesiveFrictionContact : public FrictionContact
{





  SP::SiconosVector _q_cohesion;
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CohesiveFrictionContact);

public:

  /** constructor (solver id and dimension)
      \param dimPb dimension (2D or 3D) of the friction-contact problem
      \param numericsSolverId id of the solver to be used, optional,
      default : SICONOS_FRICTION_3D_NSGS
      \rst
      see :ref:`problems_and_solvers` for details.
      \endrst
  */
  CohesiveFrictionContact(int dimPb=3, int numericsSolverId = SICONOS_FRICTION_3D_NSGS);

  /**  constructor from a pre-defined solver options set.
       \param options, the options set,
       \rst
       see :ref:`problems_and_solvers` for details.
       \endrst
  */
  CohesiveFrictionContact(int dimPb, SP::SolverOptions options);

  /** destructor
   */
  virtual ~CohesiveFrictionContact(){};


  // --- Others functions ---
  /** To compute a part of the "q+cohesion" vector of the OSNS that corresponds
   * to the cohesion (a shift the inequality constraint)
   * \param vertex, vertex (interaction) which corresponds to the considered block
   * \param pos the position of the first element of yOut to be set
  */
  void compute_q_cohesion_Block(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos);

  /** compute vector q
   *  \param time the current time
   */
  void computeq(double time);

  /** build problem coefficients (if required)
      \param time the current time
      \return true if succeeded
   */
  bool preCompute(double time);


  /** build problem coefficients (if required)
      \param time the current time
      \return true if succeeded
  */
  void postCompute();



  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

#endif // CohesiveFrictionContact_H
