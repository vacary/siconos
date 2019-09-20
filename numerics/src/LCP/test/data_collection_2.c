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

#include "lcp_test_utils.h"

char ** data_collection()
{

  int n_data_1=150;

  char ** data_collection_1 = (char **)malloc(n_data_1*sizeof(char *));
  int n_data=0;
  data_collection_1[n_data++] = "./data/lcp_enum_fails.dat";
  data_collection_1[n_data++] = "./data/lcp_CPS_1.dat";
  data_collection_1[n_data++] = "./data/lcp_CPS_2.dat";
  data_collection_1[n_data++] = "./data/lcp_CPS_3.dat";
  data_collection_1[n_data++] = "./data/lcp_CPS_4.dat";
  data_collection_1[n_data++] = "./data/lcp_CPS_4bis.dat";
  data_collection_1[n_data++] = "./data/lcp_CPS_5.dat";
  
  data_collection_1[n_data++] = "---";


  return data_collection_1;
}
