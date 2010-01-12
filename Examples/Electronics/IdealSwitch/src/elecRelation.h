#ifndef ELECRELATION_H
#define ELECRELATION_H

#include "SiconosKernel.hpp"

class elecRelation : public FirstOrderType2R
{
protected:
public:
  elecRelation();
  virtual ~elecRelation() {};


  virtual void initialize(SP::Interaction inter);


  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeh(double) ;

  /** default function to compute g
   *  \param double : current time
   */
  virtual void computeg(double) ;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJachx(double);
  virtual void computeJachlambda(double);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacgx(double);
  virtual void computeJacglambda(double);


  double source(double t);

};

TYPEDEF_SPTR(elecRelation);

#endif
