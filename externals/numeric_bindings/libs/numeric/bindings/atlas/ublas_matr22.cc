
// BLAS level 2

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/trans.hpp>
#ifdef F_USE_STD_VECTOR
#include <vector>
#include <boost/numeric/bindings/std/vector.hpp> 
#endif 
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;

using std::cout;
using std::endl; 

#ifndef F_USE_STD_VECTOR
typedef ublas::vector<double> vct_t;
typedef ublas::matrix<double, ublas::row_major> m_t;
#else
typedef ublas::vector<double, std::vector<double> > vct_t;
typedef ublas::matrix<double, ublas::column_major, std::vector<double> > m_t;
#endif 

int main() {

  cout << endl; 

  m_t m (10, 8);
  init_m (m, times_plus<double> (10, 1, 1)); 
  print_m (m, "m"); 
  cout << endl; 

  vct_t v (8);
  blas::set (1., v);
  print_v (v, "v"); 
  cout << endl; 
  vct_t vy (12); // vector size can be larger
                // than corresponding matrix size 

  // vy = m v 
  blas::gemv ( 1.0, m, v, 0.0, vy);
  print_v (vy, "vy = m v"); 
  blas::gemv (1.0, m, v, 0.0, vy);
  print_v (vy, "vy = m v"); 
  blas::gemv ( 1.0, m, v, 0.0, vy);
  print_v (vy, "vy = m v"); 
  cout << endl; 

  blas::set (1, vy); 
  print_v (vy, "vy"); 

  // v = m^T vy 
  blas::gemv ( 1.0, bindings::trans(m), vy, 0.0, v);
  print_v (v, "v = m^T vy"); 
  cout << endl; 

  // vy = 2.0 m v + 0.5 vy 
  blas::set (1., v);
  print_v (v, "v"); 
  print_v (vy, "vy"); 
  blas::gemv (2.0, m, v, 0.5, vy);
  print_v (vy, "vy = 2.0 m v + 0.5 vy"); 
  cout << endl; 
  cout << endl; 

  blas::set (1, v); 
  blas::set (0, vy); 
  print_v (v, "v"); 
  print_v (vy, "vy"); 
  cout << endl; 

  // m[2..8][1..7] 
  ublas::matrix_range<m_t> mr (m, ublas::range (2, 8), ublas::range (1, 7)); 
  print_m (mr, "mr = m[2..8][1..7]");
  cout << endl; 
  
  // vy = m[2..8][1..7] v 
  blas::gemv ( 1.0, mr, v, 0.0, vy); 
  print_v (vy, "vy = mr v"); 
  cout << endl; 

  // vy = m[2..8][1..7]^T v 
  blas::gemv ( 1.0, bindings::trans(mr), v, 0.0, vy);
  print_v (vy, "vy = mr^T v"); 
  cout << endl; 

  cout << endl; 

  // mrr = (m[2..8][1..7])[1..4][2..5]
  ublas::matrix_range<ublas::matrix_range<m_t> > 
    mrr (mr, ublas::range (1, 4), ublas::range (2, 5)); 
  print_m (mrr, "mrr = (m[2..8][1..7])[1..4][2..5]");
  cout << endl; 

  // vy = mrr v 
  blas::set (0, vy); 
  blas::gemv ( 1.0, mrr, v, 0.0, vy);
  print_v (vy, "vy = mrr v"); 
  cout << endl; 

#ifdef F_COMPILATION_FAILURE
  ublas::matrix_slice<m_t> 
    ms (m, ublas::slice (2, 1, 4), ublas::slice (1, 2, 4)); 
  print_m (ms, "ms = m[2:1:4][1:2:4]");
  cout << endl; 

  blas::set (0, vy); 
  blas::gemv ( 1.0, ms, v, 0.0, vy);
  print_v (vy, "vy = ms v"); 
  cout << endl; 
#endif 

}
