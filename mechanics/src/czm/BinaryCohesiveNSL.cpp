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

#include <algorithm>
#include <iostream>

#include "Interaction.hpp"
#include "NewtonEuler1DR.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosFwd.hpp"
#include "op3x3.h"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
//#define DEBUG_BEGIN_END_ONLY
//#define DEBUG_WHERE_MESSAGES
#include <siconos_debug.h>

// Default (private)
BinaryCohesiveNSL::BinaryCohesiveNSL() : CohesiveZoneModelNIFNSL() {}
BinaryCohesiveNSL::BinaryCohesiveNSL(unsigned int size) : CohesiveZoneModelNIFNSL(size) {}

BinaryCohesiveNSL::BinaryCohesiveNSL(double en, double et, double mu, double sigma_c,
                                     double delta_c, unsigned int size)
    : CohesiveZoneModelNIFNSL(en, et, mu, size),
      _sigma_c(sigma_c),
      _delta_c(delta_c),
      _shape_type(DOOR_SHAPE) {}

BinaryCohesiveNSL::BinaryCohesiveNSL(double en, double et, double mu, double sigma_c,
                                     double delta_c, unsigned int size,
                                     shape_type_t shape_type)
    : CohesiveZoneModelNIFNSL(en, et, mu, size),
      _sigma_c(sigma_c),
      _delta_c(delta_c),
      _shape_type(shape_type) {
  if (_shape_type == TRIANGLE_SHAPE) {
    _slope = -1.0 / _delta_c;
  }
}

BinaryCohesiveNSL::~BinaryCohesiveNSL() {}

SP::VectorOfVectors BinaryCohesiveNSL::initializeInternalVariables(Interaction& inter) {
  SP::VectorOfVectors internalVariables_sp(new VectorOfVectors());
  internalVariables_sp->resize(BinaryCohesiveNSL::INTERNAL_VARIABLE_LENGTH);

  VectorOfVectors& internalVariables = *internalVariables_sp;

  SP::Relation rel = inter.relation();
  SP::NewtonEuler1DR rel_NewtonEuler1DR(std::dynamic_pointer_cast<NewtonEuler1DR>(rel));
  if (rel_NewtonEuler1DR) {
    /* internalVariables(0) --> beta */
    /* internalVariables(1:nslawsize) --> r_cohesion */
    /* internalVariables(nslawsize+1) --> surface of the cohesive element*/
    /* Cumulative normal and tangent displacement must be also added */

    internalVariables[BinaryCohesiveNSL::R_COHESION].reset(new SiconosVector(3));
    internalVariables[BinaryCohesiveNSL::DISPLACEMENT_JUMP].reset(new SiconosVector(3));
    internalVariables[BinaryCohesiveNSL::BETA_SURFACE].reset(new SiconosVector(2));

    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])
        .setValue(0, 1.0);  // initial value of beta. This has to be fixed correctly
    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])
        .setValue(1, 1.0);  // initial value of surface. This has to be fixed correctly

    DEBUG_EXPR(std::cout << "\n The relation if of type NewtonEuler1DR" << std::endl;);

    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_1].reset(
        new SiconosVector(*rel_NewtonEuler1DR->relPc1()));
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_2].reset(
        new SiconosVector(*rel_NewtonEuler1DR->relPc2()));

    SP::SiconosVector r_nc = rel_NewtonEuler1DR->relNc();
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_NORMAL].reset(
        new SiconosVector(*r_nc));

    double Nx = r_nc->getValue(0);
    double Ny = r_nc->getValue(1);
    double Nz = r_nc->getValue(2);
    double t[6];
    double* pt = t;
    if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
      THROW_EXCEPTION(
          "NewtonEuler3DR::FC3DcomputeJachqTFromContacts. Problem in calling "
          "orthoBaseFromVector");
    pt = t;

    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1].reset(
        new SiconosVector(3));
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2].reset(
        new SiconosVector(3));

    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1]->setValue(0, pt[0]);
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1]->setValue(1, pt[1]);
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1]->setValue(2, pt[2]);
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2]->setValue(0, pt[3]);
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2]->setValue(1, pt[4]);
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2]->setValue(2, pt[5]);

    SP::SiconosVector pc1 = rel_NewtonEuler1DR->pc1();
    SP::SiconosVector pc2 = rel_NewtonEuler1DR->pc2();
    SP::SiconosVector nc = rel_NewtonEuler1DR->nc();

    internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_1].reset(new SiconosVector(*pc1));
    internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_2].reset(new SiconosVector(*pc2));
    internalVariables[BinaryCohesiveNSL::NORMAL].reset(new SiconosVector(*nc));

    internalVariables[BinaryCohesiveNSL::TANGENT_1].reset(new SiconosVector(3));
    internalVariables[BinaryCohesiveNSL::TANGENT_2].reset(new SiconosVector(3));

    Nx = nc->getValue(0);
    Ny = nc->getValue(1);
    Nz = nc->getValue(2);
    if (orthoBaseFromVector(&Nx, &Ny, &Nz, pt, pt + 1, pt + 2, pt + 3, pt + 4, pt + 5))
      THROW_EXCEPTION(
          "NewtonEuler3DR::FC3DcomputeJachqTFromContacts. Problem in calling "
          "orthoBaseFromVector");
    pt = t;

    internalVariables[BinaryCohesiveNSL::TANGENT_1]->setValue(0, pt[0]);
    internalVariables[BinaryCohesiveNSL::TANGENT_1]->setValue(1, pt[1]);
    internalVariables[BinaryCohesiveNSL::TANGENT_1]->setValue(2, pt[2]);
    internalVariables[BinaryCohesiveNSL::TANGENT_2]->setValue(0, pt[3]);
    internalVariables[BinaryCohesiveNSL::TANGENT_2]->setValue(1, pt[4]);
    internalVariables[BinaryCohesiveNSL::TANGENT_2]->setValue(2, pt[5]);

    SiconosVector displacement_jump(3);
    displacement_jump = *pc2 - *pc1;
    DEBUG_EXPR(displacement_jump.display(););
    internalVariables[BinaryCohesiveNSL::INITIAL_DISPLACEMENT_JUMP].reset(
        new SiconosVector(displacement_jump));
    DEBUG_EXPR(for (auto v
                    : internalVariables) {
      if (v) v->display();
    };);

  } else {
    /* internalVariables(0) --> beta */
    /* internalVariables(1:nslawsize) --> r_cohesion */
    /* internalVariables(nslawsize+1) --> surface of the cohesive element*/
    /* Cumulative normal and tangent displacement must be also added */
    internalVariables[BinaryCohesiveNSL::R_COHESION].reset(new SiconosVector(3));

    internalVariables[BinaryCohesiveNSL::BETA_SURFACE].reset(new SiconosVector(2));

    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])
        .setValue(0, 1.0);  // initial value of beta. This has to be fixed correctly
    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])
        .setValue(1, 1.0);  // initial value of surface. This has to be fixed correctly
  }
  // getchar();
  return internalVariables_sp;
}
void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter) {
  DEBUG_BEGIN("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");

  VectorOfVectors& internalVariables = *inter.internalVariables();
  VectorOfVectors& internalVariables_k = *inter.internalVariables_k();

  // double * internalVariablesArray =  inter.internalVariables()->getArray();
  // double * beta_k = &(inter.internalVariables_k()->getArray()[0]);
  // double  surface = inter.internalVariables_k()->getArray()[_size+1];

  double* beta = &(internalVariables[BinaryCohesiveNSL::BETA_SURFACE]->getArray()[0]);
  double* surface = &(internalVariables[BinaryCohesiveNSL::BETA_SURFACE]->getArray()[1]);

  double* beta_k = &(internalVariables_k[BinaryCohesiveNSL::BETA_SURFACE]->getArray()[0]);

  double u_N = 0.0;
  double u_T = 0.0;
  double u_S = 0.0;

  if (*beta_k > 0.0) {
    /*  compute displacement */

    double delta = 0;
    SP::Relation rel = inter.relation();
    SP::NewtonEuler1DR rel_NewtonEuler1DR(std::dynamic_pointer_cast<NewtonEuler1DR>(rel));
    if (rel_NewtonEuler1DR) {
      SiconosVector& r_pc1_0 =
          *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_1];
      SiconosVector& r_pc2_0 =
          *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_2];
      SiconosVector& r_nc_0 = *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_NORMAL];
      SiconosVector& r_t1_0 =
          *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1];
      SiconosVector& r_t2_0 =
          *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2];
      SiconosVector& pos_0 = *internalVariables[BinaryCohesiveNSL::DISPLACEMENT_JUMP];

      // Compute the current cohesive point in absolute frame.
      VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
      BlockVector& q = *DSlink[NewtonEulerR::q0];
      SiconosVector& pc1_0 = *internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_1];
      SiconosVector& pc2_0 = *internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_2];
      SiconosVector& nc_0 = *internalVariables[BinaryCohesiveNSL::NORMAL];
      SiconosVector& t1_0 = *internalVariables[BinaryCohesiveNSL::TANGENT_1];
      SiconosVector& t2_0 = *internalVariables[BinaryCohesiveNSL::TANGENT_2];

      rel_NewtonEuler1DR->computeContactPointsFromRelativeContactPoints(
          q, r_pc1_0, r_pc2_0, r_nc_0, r_t1_0, r_t2_0, pc1_0, pc2_0, nc_0, t1_0, t2_0);

      DEBUG_EXPR(std::cout << "pc1_0 is: " << pc1_0 << std::endl;
                 std::cout << "pc2_0 is: " << pc2_0 << std::endl;
                 std::cout << "nc_0 is: " << nc_0 << std::endl;);

      SiconosVector pos(3);
      pos = pc2_0 - pc1_0;

      // WARNING: curious sign. may depend on which is the ds2?
      pos = pc1_0 - pc2_0;

      DEBUG_EXPR(std::cout << "pos :" << pos << std::endl;);
      SiconosVector u(3);
      u = pos - pos_0;

      u_N = inner_prod(u, nc_0);
      u_T = inner_prod(u, t1_0);
      u_S = inner_prod(u, t2_0);

      DEBUG_EXPR(std::cout << "displacement jump u :" << u << std::endl;);
      delta = u.norm2();
      DEBUG_EXPR(std::cout << "delta :" << delta << std::endl;);

    } else {
      if (_size == 1) {
        delta = (*(inter.y(0)))(
            0);  // this rule has to be improved following the model of Tveergard.
      } else {
        delta = inter.y(0)->norm2();
        SiconosVector& y = *(inter.y(0));

        //	      y.display();

        u_T = y(1);
        u_S = y(2);
      }
    }

    /*  update beta */

    // std::cout << this << std::endl;
    DEBUG_PRINTF("beta = %e\n", *beta);
    DEBUG_PRINTF("beta_k = %e\n", *beta_k);
    DEBUG_PRINTF("delta = %e\n", delta);

    if (_shape_type == DOOR_SHAPE) {
      if ((delta > _delta_c)) {
        DEBUG_PRINT("the interface is broken\n");
        *beta = 0.0;
      } else if ((delta <= _delta_c) and (*beta_k == 1.0)) {
        DEBUG_PRINT("the interface is sane\n");
        *beta = 1.0;
      }
    } else if (_shape_type == TRIANGLE_SHAPE) {
      *beta = std::min(*beta_k, 1.0 + _slope * delta);
      *beta = std::max(0., *beta);
    }
  } else {
    *beta = *beta_k;
  }

  DEBUG_PRINTF("beta = %e\n", *beta);
  double* r_cohesion = internalVariables[BinaryCohesiveNSL::R_COHESION]->getArray();

  /* compute _r_cohesion */
  for (int k = 1; k < _size; k++) {
    r_cohesion[k] = 0.0;
  }
  r_cohesion[0] = -*beta * _sigma_c * *surface;

  DEBUG_PRINTF("normal  cohesion force %4.2e\n", r_cohesion[0]);
  /* compute explit tangential part of _r_cohesion */

  // r_cohesion[1]= - *beta * _sigma_c * surface; not possible with an extrinsic cohesive law
  if (_size > 2) {
    double norm_u_T = sqrt(u_T * u_T + u_S * u_S);

    if (norm_u_T > 0.0) {
      double d_T1 = u_T / norm_u_T;
      double d_T2 = u_S / norm_u_T;

      r_cohesion[1] = -*beta * _sigma_c * *surface * d_T1;
      r_cohesion[2] = -*beta * _sigma_c * *surface * d_T2;

      DEBUG_PRINTF("tangential cohesion force %4.2e\t %4.2e\n", r_cohesion[1], r_cohesion[2]);
    } else {
      r_cohesion[1] = 0.0;
      r_cohesion[2] = 0.0;
    }
  }

  DEBUG_EXPR(std::cout << "r_cohesion is: "
                       << *internalVariables[BinaryCohesiveNSL::R_COHESION] << std::endl;);

  // getchar();
  assert(*beta <= *beta_k);

  DEBUG_END("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");
}
bool BinaryCohesiveNSL::isActiveAtLevel(Interaction& inter, unsigned int level) {
  // double * beta = &(inter.internalVariables()->getArray()[0]);
  // if (level <=1)
  // {
  //   if (*beta > 0.0)
  //   {
  //     return true; // when the interface is cohesive, we force the activation of the
  //     constraint at the veloicity level
  //   }
  //   else
  //   {
  //     return false;
  //   }
  // }
  // else
  //   THROW_EXCEPTION("BinaryCohesiveNSL::isActiveAtLevel(unsigned int level): level should be
  //   less than 1");
  return false;
}

double* BinaryCohesiveNSL::r_cohesion(Interaction& inter) const {
  VectorOfVectors& internalVariables = *inter.internalVariables();
  double* r_cohesion = internalVariables[BinaryCohesiveNSL::R_COHESION]->getArray();
  return r_cohesion;
};

double BinaryCohesiveNSL::beta(Interaction& inter) const {
  VectorOfVectors& internalVariables = *inter.internalVariables();
  double beta = internalVariables[BinaryCohesiveNSL::BETA_SURFACE]->getValue(0);
  return beta;
};

void BinaryCohesiveNSL::display() const {
  CohesiveZoneModelNIFNSL::display();
  std::cout << "=== BinaryCohesiveNSL data display ===============================" << this
            << std::endl;
  std::cout << " cohesive resistance to traction: " << _sigma_c << std::endl;
  std::cout << " critical displacement: " << _delta_c << std::endl;
  std::cout << "=================================================================="
            << std::endl;
}
