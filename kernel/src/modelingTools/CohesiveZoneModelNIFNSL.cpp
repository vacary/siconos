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
#include "CohesiveZoneModelNIFNSL.hpp"

#include <iostream>
// Default (private)
CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL():
  NewtonImpactFrictionNSL()
{
}
CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL(unsigned int size):
  NewtonImpactFrictionNSL(size , 0.0, 0.0, 0.0)
{}

CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL(double en, double et, double mu, unsigned int size):
  NewtonImpactFrictionNSL(en, et, mu, size)
{
  _r_cohesion = new double[size];
  for (int k=0 ; k < size ; k++) _r_cohesion[k]=0.0;
}

CohesiveZoneModelNIFNSL::~CohesiveZoneModelNIFNSL()
{}


void CohesiveZoneModelNIFNSL::display() const
{
  NewtonImpactFrictionNSL::display();
  std::cout << "=== CohesiveZoneModelNIFNSL data display ===============================" << std::endl;
  std::cout << " normal current cohesive reaction: " << _r_cohesion[0] <<std::endl;
  std::cout << "==================================================================" <<std::endl;
}
