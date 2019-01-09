/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <errno.h>
char *** test_collection(int, char **);

char *** test_collection(int n_data_1, char ** data_collection)
{
  int n_test=150;
  int n_entry = 50;
  char *** test_nsgs = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test_nsgs[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  int e=0;
  int d=8;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "---";
  n++;
  
  e=0;
  d=9;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "---";
  n++;
  
  test_nsgs[n][0] ="---";
  return test_nsgs;
  
}
