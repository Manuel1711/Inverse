#include <boost/math/quadrature/gauss_kronrod.hpp> 
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "gmpfrxx.cpp" 
#include <eigen3/Eigen/Dense>
#include <vector>
    
#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

using namespace std;
using namespace Eigen;
//using T=mpfr_class;

//Import of the boost multiprecision libraries
namespace bm=boost::multiprecision;
namespace bq=boost::math::quadrature;

//Define of a arbitrary precise variable Real
using Real=
  bm::number<bm::mpfr_float_backend<1024>>;
const Real inf=
  std::numeric_limits<Real>::infinity();

//Using the eigen library, define arbitrary precise
//matrix and vectors
using PrecMatr= Matrix<Real,Dynamic,Dynamic>;
using PrecVec = Matrix<Real,Dynamic, 1>;


//Functions for printing as a string the Real variable
string conv(const Real& in)
{
  ostringstream os;
  os.precision(mpf_get_default_prec()/4);
  
  os<<in;

  return os.str();
}
 
Real conv(const string& in)
{
  istringstream is(in);
  
  Real out;
  
  is>>out;
  
  return out;
}

