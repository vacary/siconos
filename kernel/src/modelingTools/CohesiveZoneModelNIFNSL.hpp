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
/*! \file BinaryCohesiveNSL.hpp
  Binary cohesive model with impat and friction

  The model implements a naive binary (broken or sane) cohesive zone model.
  When the interface is sane, the limt traction threshols is constant and equal to sigma_c.
  The interface is broken if the nomral displacement id larger than delta_c
*/

#ifndef CohesiveZoneModelNIFNSL_H
#define CohesiveZoneModelNIFNSL_H

#include "NewtonImpactFrictionNSL.hpp"

/** CohesiveZoneModelNIFNSL
 * Base class for cohesive zone model based on  NewtonImpactFrictionNSL
 */
class CohesiveZoneModelNIFNSL : public NewtonImpactFrictionNSL
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CohesiveZoneModelNIFNSL);


protected:

  /** default constructor
   */
  CohesiveZoneModelNIFNSL();

public:

  /** basic constructor
   *  \param size size of the ns law
   */
  explicit CohesiveZoneModelNIFNSL(unsigned int size);

  /** constructor with the value of the CohesiveZoneModelNIFNSL attributes
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param sigma_c double : cohesive resistance to traction 
   *  \param delta_c double : critical displacement 
   *  \param size unsigned int: size of the ns law
   */
  CohesiveZoneModelNIFNSL(double en, double et, double mu, unsigned int size);

  /** Destructor */
  ~CohesiveZoneModelNIFNSL();

  /** getter of r_cohesion
   * \return the value of r_cohesion
   */
  virtual double * r_cohesion(Interaction& inter) const =0;
  
  // OTHER FUNCTIONS

  /** print the data to the screen
   */
  void display() const;

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();

};
DEFINE_SPTR(CohesiveZoneModelNIFNSL)
#endif // CohesiveZoneModelNIFNSL_H
