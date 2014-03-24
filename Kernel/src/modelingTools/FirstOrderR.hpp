/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file FirstOrderR.hpp
\brief General interface for relations.
 */

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.hpp"
#include "Interaction.hpp"

/** Pointer to function for plug-in for operators related to output and its gradients.*/
typedef void (*OutPtr)(unsigned int, double*, double, unsigned int, double*, double*, unsigned int, double*);

/** Pointer to function for plug-in for operators related to input and its gradients.*/
typedef void (*InPtr)(unsigned int, double*, double, unsigned int, double*, unsigned int, double*);

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(X,t,\lambda,Z)\\
 * R &=& g(X,t,\lambda,Z)
 * \f}
 *  X, Z, R corresponds to DynamicalSystem variables.
 *  If DS1 and DS2 are involved in the linked Interaction, then X =[x1 x2], Z=[z1 z2] ...
 *
 *  \f$ y \ and \ \lambda \f$ are specific variables of the interaction (see this class for more details).
 *  h and g are plugged on external functions, via plug-in mechanism (see SiconosSharedLibrary).
 *
 * h <=> output
 *
 * g <=> input
 *
 * Operators (and their corresponding plug-in):
- h: saved in Interaction as y (plug-in: output[0])
- \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
- \f$ \nabla_\lambda h \f$: jacobianH[1] ( output[2] )
- g: saved in DS as r ( input[0])
- \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )


Note: we use a vector for jacobianG while there is only one jacobian. Just for future changes and to allow easy new implementations if some other
variables are required in g.

 *
 */
class FirstOrderR : public Relation
{
public:

  enum DataNames {free, z, x, xq, r, deltax, g_alpha, residu_r, ds_xp, sizeDataNames};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderR);

  /** basic constructor
  *  \param the type of the relation
  */
  FirstOrderR(RELATION::SUBTYPES newType): Relation(RELATION::FirstOrder, newType) {}

  SP::SiconosMatrix _jachx;
  SP::SiconosMatrix _jachz;
  SP::SiconosMatrix _jacglambda;

public:

  /** destructor
  */
  virtual ~FirstOrderR() {};

  /** get a pointer on matrix Jacg[index]
  *  \return a pointer on a SiconosMatrix
  */
  virtual SP::SiconosMatrix jacglambda() const
  {
    return _jacglambda;
  }

  /** set Jacg[index] to pointer newPtr (pointer link)
  *  \param SP::SiconosMatrix  newPtr
  *  \param unsigned int: index position in Jacg vector
  */
  inline void setJacglambdaPtr(SP::SiconosMatrix newPtr)
  {
    _jacglambda = newPtr ;
  }

  /** initialize the relation (check sizes, memory allocation ...)
  \param SP to Interaction: the interaction that owns this relation
  */
  virtual void initialize(Interaction& inter);

  /** default function to compute h
  *  \param double : current time
  */
  virtual void computeh(double time, Interaction& inter) = 0;

  /** default function to compute g
  *  \param double : current time
  */
  virtual void computeg(double time, Interaction& inter) = 0;

  /** default function to compute jacobianH
  *  \param double : current time
  *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
  */
  virtual void computeJachx(double time, Interaction& inter);
  virtual void computeJachz(double time, Interaction& inter);
  virtual void computeJachlambda(double time, Interaction& inter);
  virtual void computeJach(double time, Interaction& inter)
  {
    computeJachx(time, inter);
    computeJachz(time, inter);
    computeJachlambda(time, inter);
  }


  /** default function to compute jacobianG according to lambda
  *  \param double : current time
  *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
  */
  virtual void computeJacglambda(double time, Interaction& inter);

  virtual void computeJacg(double time, Interaction& inter)
  {
    computeJacglambda(time, inter);
  }

  SP::SiconosMatrix jachx() const
  {
    return _jachx;
  }

  /**
  * return a SP on the C matrix.
  * The matrix C in the linear case, else it returns Jacobian of the output with respect to x.
  *
  */
  inline SP::SiconosMatrix C() const
  {
    return _jachx;
  }
  /**
  * return a SP on the F matrix.
  * The matrix F in the linear case, else it returns Jacobian of the output with respect to z.
  *
  */
  inline SP::SiconosMatrix F() const
  {
    return _jachz;
  }
  /**
  * return a SP on the D matrix.
  * The matrix D in the linear case, else it returns Jacobian of the output with respect to lambda.
  */
  inline SP::SiconosMatrix D() const
  {
    return _jachlambda;
  }
  /**
  * return a SP on the B matrix.
  * The matrix B in the linear case, else it returns Jacobian of the input with respect to lambda.
  */
  inline SP::SiconosMatrix B() const
  {
    return _jacglambda;
  }

};
TYPEDEF_SPTR(FirstOrderR)

#endif
