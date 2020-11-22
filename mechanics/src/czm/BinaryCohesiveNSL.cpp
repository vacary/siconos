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

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
//#define DEBUG_BEGIN_END_ONLY
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

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
  _sigma_c(sigma_c), _surface(1.0),
  _delta_c(delta_c), _beta(1.0), _beta_k(1.0)
{}



BinaryCohesiveNSL::~BinaryCohesiveNSL()
{}

SP::SiconosVector BinaryCohesiveNSL::initializeInternalVariables(Interaction& inter)
{
  /* internalVariables(0) --> beta */
  /* internalVariables(1:nslawsize) --> r_cohesion */
  SP::SiconosVector internalVariables(new SiconosVector(1+_size));
  
  return internalVariables;
  
}
void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter) 
{
  DEBUG_BEGIN("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");
  /*  update beta */
  
  double normal_gap = (*(inter.y(0)))(0); // this rule has to be improved following the model of Tveergard.
  // std::cout << this << std::endl;
  DEBUG_PRINTF("beta_k = %e\n", _beta_k);
  DEBUG_PRINTF("normal_gap = %e\n", normal_gap);
  
  if ((normal_gap > _delta_c))
  {
    DEBUG_PRINT("the interface is broken\n");
    _beta=0.0;
  }
  else if ((normal_gap <= _delta_c) and (_beta_k == 1.0))
  {
    DEBUG_PRINT("the interface is sane\n");
    _beta=1.0;
  }
  
  DEBUG_PRINTF("beta = %e\n", _beta);

  double * r_cohesion = &(inter.internalVariables()->getArray()[1]);
  /* compute _r_cohesion */
  for (int k =1; k < _size; k++)
  {
    r_cohesion[k]= 0.0;
  }
  r_cohesion[0]= - _beta* _sigma_c * _surface;

  assert(_beta <= _beta_k);
  
  DEBUG_END("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");
}
bool BinaryCohesiveNSL::isActiveAtLevel(Interaction& inter, unsigned int level)
{

  if (level <=1)
  {
    if (_beta > 0.0)
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


void BinaryCohesiveNSL::display() const
{
  CohesiveZoneModelNIFNSL::display();
  std::cout << "=== BinaryCohesiveNSL data display ===============================" << this << std::endl;
  std::cout << " cohesive resistance to traction: " << _sigma_c <<std::endl;
  std::cout << " cohesive surface: " << _surface <<std::endl;
  std::cout << " critical displacement: " << _delta_c <<std::endl;
  std::cout << " cohesion state: " << _beta <<std::endl;
  std::cout << " old cohesion state: " << _beta_k <<std::endl;
  std::cout << "==================================================================" <<std::endl;
}
