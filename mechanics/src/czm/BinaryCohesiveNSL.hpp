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

#ifndef BINARYCOHESIVENSLAW_H
#define BINARYCOHESIVENSLAW_H

#include "CohesiveZoneModelNIFNSL.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosFwd.hpp"

/** Newton Impact-Friction Non Smooth Law
 *
 */
class BinaryCohesiveNSL : public CohesiveZoneModelNIFNSL
{

public:
  typedef enum {DOOR_SHAPE,TRIANGLE_SHAPE} shape_type_t;
  enum BinaryCohesiveNSL_internalVariables {
    BETA_SURFACE,
    R_COHESION,
    DISPLACEMENT_JUMP,
    COHESIVE_POINT_1,
    COHESIVE_POINT_2,
    NORMAL,
    TANGENT_1,
    TANGENT_2,
    INITIAL_DISPLACEMENT_JUMP,
    INITIAL_RELATIVE_COHESIVE_POINT_1,
    INITIAL_RELATIVE_COHESIVE_POINT_2,
    INITIAL_RELATIVE_NORMAL,
    INITIAL_RELATIVE_TANGENT_1,
    INITIAL_RELATIVE_TANGENT_2,
    INTERNAL_VARIABLE_LENGTH
  };

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BinaryCohesiveNSL);

  /** cohesive resistance to traction */
  double _sigma_c;

  /** critical displacement */
  double _delta_c;

  /** slope for triangle law */
  double _slope;

  /** shape of the cohesive law */
  shape_type_t _shape_type;

  /** fallback law when the interface is broken */
  SP::NonSmoothLaw _nslaw_broken;


protected:

  /** default constructor
   */
  BinaryCohesiveNSL();

public:

  /** basic constructor
   *  \param size size of the ns law
   */
  explicit BinaryCohesiveNSL(unsigned int size);

  /** constructor with the value of the BinaryCohesiveNSL attributes
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param sigma_c double : cohesive resistance to traction
   *  \param delta_c double : critical displacement
   *  \param size unsigned int: size of the ns law
   */
  BinaryCohesiveNSL(double en, double et, double mu, double sigma_c, double delta_c, unsigned int size);

  BinaryCohesiveNSL(double en, double et, double mu,
		    double sigma_c, double delta_c,
		    unsigned int size,
		    shape_type_t shape_type);

  /** Destructor */
  ~BinaryCohesiveNSL();

  /** getter of sigma_c
   * \return the value of sigma_c
   */
  inline double sigma_c() const
  {
    return _sigma_c;
  };

  /** setter of sigma_c
   * \param newVal a double to set sigma_c
   */
  inline void setSigma_C(double newVal)
  {
    _sigma_c = newVal;
  };

  /** getter of delta_c
   * \return the value of delta_c
   */
  inline double delta_c() const
  {
    return _delta_c;
  };

  /** setter of delta_c
   * \param newVal a double to set delta_c
   */
  inline void setDelta_C(double newVal)
  {
    _delta_c = newVal;
  };


  // OTHER FUNCTIONS
  SP::VectorOfVectors initializeInternalVariables(Interaction &) override;

  void updateInternalVariables(Interaction & inter) override;

  /** Ask if the Nslaw is active at a given level
  */
  virtual bool isActiveAtLevel(Interaction& inter,  unsigned int level) override;

  double * r_cohesion(Interaction& inter) const override;

  /** getter of beta
   * \return the value of beta
   */
  double beta(Interaction& inter) const;

  /** print the data to the screen
   */
  void display() const override;

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();

};

#endif // BINARYCOHESIVENSLAW_H
