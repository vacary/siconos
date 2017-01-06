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

#ifndef TLSDEF_H
#define TLSDEF_H

/*! \file tlsdef.h
 *  \brief definition of thread local variable
 */

#if defined(__GNUC__)
#define DESTRUCTOR_ATTR __attribute__ ((destructor))
#else
#define DESTRUCTOR_ATTR 
#endif

#ifndef __cplusplus

  #if __STDC_VERSION__ >= 201112L
    #include <threads.h>
    #define tlsvar thread_local
  #else

    #if defined(__GNUC__)
      #define tlsvar __thread 
    #else
      #error "Don't know how to create a thread-local variable"
    #endif
  #endif


#else

  #if SICONOS_CXXVERSION >= 201103L
    #define tlsvar thread_local
  #else
    #if defined(__GNUC__)
      #define tlsvar __thread
    #elif defined(_MSC_VER)
      #define tlsvar __declspec(thread)
    #else
      #error "Don't know how to create a thread-local variable"
    #endif
  #endif

#endif

#endif