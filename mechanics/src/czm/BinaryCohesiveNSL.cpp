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
#include "BinaryCohesiveNSL.hpp"
#include "Interaction.hpp"

#include <iostream>
#include <algorithm>
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
//#define DEBUG_BEGIN_END_ONLY
//#define DEBUG_WHERE_MESSAGES
#include <siconos_debug.h>

// Default (private)
BinaryCohesiveNSL::BinaryCohesiveNSL():
  CohesiveZoneModelNIFNSL()
{}
BinaryCohesiveNSL::BinaryCohesiveNSL(unsigned int size):
  CohesiveZoneModelNIFNSL(size)
{}


BinaryCohesiveNSL::BinaryCohesiveNSL(double en, double et, double mu,
                                     double sigma_c, double delta_c,
                                     unsigned int size):
  CohesiveZoneModelNIFNSL(en, et, mu, size),
  _sigma_c(sigma_c),
  _delta_c(delta_c),
  _shape_type(DOOR_SHAPE)
{}

BinaryCohesiveNSL::BinaryCohesiveNSL(double en, double et, double mu,
                                     double sigma_c, double delta_c,
                                     unsigned int size,
				     shape_type_t shape_type):
  CohesiveZoneModelNIFNSL(en, et, mu, size),
  _sigma_c(sigma_c),
  _delta_c(delta_c),
  _shape_type(shape_type)
{
	if (_shape_type == TRIANGLE_SHAPE){
		_slope = - 1.0/_delta_c;
	}
}



BinaryCohesiveNSL::~BinaryCohesiveNSL()
{}

SP::SiconosVector BinaryCohesiveNSL::initializeInternalVariables(Interaction& inter)
{
  /* internalVariables(0) --> beta */
  /* internalVariables(1:nslawsize) --> r_cohesion */
  /* internalVariables(nslawsize+1) --> surface od the cohesive element*/	
  /* Cumulative normal and tangent displacement must be also added */
  SP::SiconosVector internalVariables(new SiconosVector(1+_size+1));
  internalVariables->setValue(0,1.0); // initial value of beta. This has to be fixed correctly
  internalVariables->setValue(1+_size,1.0); // initial value of surface. This has to be fixed correctly
  return internalVariables;

}
void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)
{
  DEBUG_BEGIN("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");
  /*  update beta */

  double delta = (*(inter.y(0)))(0); // this rule has to be improved following the model of Tveergard.

  double * beta = &(inter.internalVariables()->getArray()[0]);
  double * beta_k = &(inter.internalVariables_k()->getArray()[0]);
  double  surface = inter.internalVariables_k()->getArray()[_size+1];

  // std::cout << this << std::endl;
  DEBUG_PRINTF("beta = %e\n", *beta);
  DEBUG_PRINTF("beta_k = %e\n", *beta_k);
  DEBUG_PRINTF("delta = %e\n", delta);

  if (_shape_type == DOOR_SHAPE)
    {
      if ((delta > _delta_c))
	{
	  DEBUG_PRINT("the interface is broken\n");
	  *beta=0.0;
	}
      else if ((delta <= _delta_c) and (*beta_k == 1.0))
	{
	  DEBUG_PRINT("the interface is sane\n");
	  *beta=1.0;
	}
    }
  else if (_shape_type == TRIANGLE_SHAPE)
  {	  
	  *beta = std::min(*beta_k, 1.0 + _slope * delta);
	  *beta = std::max(0., *beta);	  
  }
  DEBUG_PRINTF("beta = %e\n", *beta);

  double * r_cohesion = &(inter.internalVariables()->getArray()[1]);
  /* compute _r_cohesion */
  for (int k =1; k < _size; k++)
  {
    r_cohesion[k]= 0.0;
  }
  r_cohesion[0]= - *beta * _sigma_c * surface;
  r_cohesion[1]= 0.0;
  //r_cohesion[1]= - *beta * _sigma_c * surface; not possible with an extrinsic cohesive law
  if (_size > 2)
  {
    r_cohesion[2]= r_cohesion[1];
  }
  
  assert( *beta <= *beta_k);

  DEBUG_END("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");
}
bool BinaryCohesiveNSL::isActiveAtLevel(Interaction& inter, unsigned int level)
{

  double * beta = &(inter.internalVariables()->getArray()[0]);
  if (level <=1)
  {
    if (*beta > 0.0)
    {
      return true; // when the interface is cohesive, we force the activation of the constraint at the veloicity level
    }
    else
    {
      return false;
    }
  }
  else
    THROW_EXCEPTION("BinaryCohesiveNSL::isActiveAtLevel(unsigned int level): level should be less than 1");

}

double * BinaryCohesiveNSL::r_cohesion(Interaction& inter) const
{
  double * r_cohesion = &(inter.internalVariables()->getArray()[1]);
  return r_cohesion;
};

double  BinaryCohesiveNSL::beta(Interaction& inter) const
{
  double beta = (inter.internalVariables()->getArray()[0]);
  return beta;
};


void BinaryCohesiveNSL::display() const
{
  CohesiveZoneModelNIFNSL::display();
  std::cout << "=== BinaryCohesiveNSL data display ===============================" << this << std::endl;
  std::cout << " cohesive resistance to traction: " << _sigma_c <<std::endl;
  std::cout << " critical displacement: " << _delta_c <<std::endl;
  std::cout << "==================================================================" <<std::endl;
}
