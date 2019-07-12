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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "sanitizer.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

#include "NumericsVector.h"
void NV_display(double * m, int nRow)
{
  int lin;
  printf("vector of size\t%d\t =\n[", nRow);
  if (nRow == 0)
  {
    printf("]\n");
  }
  for (lin = 0; lin < nRow; lin++)
  {
    printf(" %.15e", m[lin]);
    if (lin != nRow - 1)
      printf(", ");
    else
      printf("]\n");
  }

}
void NV_write_in_file_python(double * m,  int nRow, FILE* file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NV_write_in_file_python  failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "size = %d; \n", nRow);
  fprintf(file, "data= [");
  for (int i = 0; i < nRow; i++)
  {
    fprintf(file, "%32.24e,\t ", m[i]);
  }
  fprintf(file, "]");
}

bool NV_equal(double * x, double * y, int n, double tol)
{
  for (int i =0; i< n ; i++)
  {
    if (fabs(x[i] - y[i]) >= tol)
    {
      DEBUG_PRINTF("error %i = %e\n",i, fabs(x[i]) - y[i]);
      return false;
    } 
  }
  return true;
}

void NV_insert(double * x, const unsigned int xSize,
               const double * const y, const unsigned int ySize,
               unsigned int i)
{
    if (xSize < ySize)
    {
        fprintf(stderr, "NV_insert ::  the vector to be inserted is greater than the given vector: size_x < size_y - %d < %d\n", xSize, ySize);
        exit(EXIT_FAILURE);
    }
    if (i + ySize > xSize)
    {
        fprintf(stderr, "NV_insert ::  the vector to be inserted is too big for insertion from position %d\n", i);
        exit(EXIT_FAILURE);
    }
    for (unsigned int j = i; j < i + ySize; ++j)
        x[j] = y[j - i];
}

double* NV_power2(const double * const vec, const unsigned int vecSize)
{
    double * out = (double*)malloc(vecSize * sizeof(double));
    for (unsigned int i = 0; i < vecSize; ++i)
        out[i] = vec[i] * vec[i];
    return out;
}


double NV_reduce(const double * const vec, const unsigned int vecSize)
{
    register double sum = 0.0;
    for (unsigned int i = 0; i < vecSize; ++i)
        sum += vec[i];
    return sum;
}

double* NV_prod(const double * const x, const double * const y, const unsigned int vecSize)
{
    double * out = (double*)malloc(vecSize * sizeof(double));
    for (unsigned int i = 0; i < vecSize; ++i)
        out[i] = x[i] * y[i];
    return out;
}

double NV_min(const double * const vec, const unsigned int vecSize)
{
    double min_elem = DBL_MAX;
    for (unsigned int i = 0; i < vecSize; ++i)
        if (vec[i] < min_elem)
            min_elem = vec[i];
    return min_elem;
}
