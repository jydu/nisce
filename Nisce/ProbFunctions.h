/*
 * =====================================================================================
 *
 *       Filename:  ProbFunctions.h
 *
 * =====================================================================================
 */

#ifndef _PROB_FUNCTIONS_H_
#define _PROB_FUNCTIONS_H_

// From the STL:
#include <iterator>

// From Boost:
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace std;

using d1_array_ui = vector<unsigned int>;
using d2_array_ui = vector<d1_array_ui>;
using d1_array_dbl = vector<double>;
using d2_array_dbl = vector<d1_array_dbl>;
using d3_array_dbl = vector<d2_array_dbl>;
using d4_array_dbl = vector<d3_array_dbl>;
using d1_array_str = vector<string>;

// LogSumExp function, https://stackoverflow.com/questions/45943750/calculating-log-sum-exp-function-in-c
template<typename Iter>
typename iterator_traits<Iter>::value_type logSumExp(Iter begin, Iter end)
{
  using VT = typename iterator_traits<Iter>::value_type;
  if (begin == end)
    return VT{}
  ;
  auto max_elem = *max_element(begin, end);
  auto sum = accumulate(begin, end, VT{},
                        [max_elem](VT a, VT b) {
    return a + exp(b - max_elem);
  });
  return max_elem + log(sum);
}

double dpoisson(unsigned int x, double lambda)
{
  if (lambda == 0)
    return x == 0 ? 1.0 : 0.0;
  if (x < 0)
    return 0.;
  if (x <= lambda * DBL_MIN)
    return exp(-lambda);
  double p = exp(-lambda + x * log(lambda) - lgamma(x + 1));
  return max(p, std::numeric_limits< double >::min());
}

template<class T>
class NumberResult
{
public:
  T value;
  T score;

  NumberResult(T v, T s):
    value(v),
    score(s)
  {}
};

#endif //_PROB_FUNCTIONS_H_

