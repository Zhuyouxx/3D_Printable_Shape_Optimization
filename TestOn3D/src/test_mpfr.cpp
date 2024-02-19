//#include <stdlib.h>
//#include <math.h>
//#include <string>
//#include <boost/multiprecision/mpfr.hpp>
//
//
//int main() {
//    boost::multiprecision::mpfr_float::default_precision(1000);
//    typedef boost::multiprecision::mpfr_float mpfr_float;
//    // Define two high-precision numbers
//    mpfr_float a = 123456789.123456789123456789;
//    mpfr_float b = 987654321.987654321987654321;
//    // Perform arithmetic operations
//    mpfr_float sum = a + b;
//    mpfr_float product = a * b;
//
//    // Set the precision for output to 50 digits
//    std::cout.precision(100);
//
//    // Output the results
//    std::cout << "Sum of a and b: " << sum << std::endl;
//    std::cout << "Product of a and b: " << product << std::endl;
//
//
//
//    double a1 = 123456789.123456789123456789;
//    double b1 = 987654321.987654321987654321;
//    // Perform arithmetic operations
//    double sum1 = a1 + b1;
//    double product1 = a1 * b1;
//
//    // Set the precision for output to 50 digits
//    std::cout.precision(50);
//
//    // Output the results
//    std::cout << "Sum of a and b: " << sum1 << std::endl;
//    std::cout << "Product of a and b: " << product1 << std::endl;
//
//    return 0;
//}

#include <iostream>
#include <boost/multiprecision/mpfr.hpp>
//#include <boost/math/special_functions/gamma.hpp>

# include <cppad/cppad.hpp> // the CppAD package
# include <cppad/example/cppad_eigen.hpp>
# include <cppad/example/atomic_two/eigen_mat_inv.hpp>
# include <cppad/example/atomic_two/eigen_mat_mul.hpp>

// Define the namespace for convenience
namespace mp = boost::multiprecision;
int main() {
    using namespace boost::multiprecision;
    using ADmpfr_float = CppAD::AD<mpfr_float>;
    // Operations at variable precision and no numeric_limits support:
    using namespace std;
    mpfr_float a = 2;
    mpfr_float::default_precision(1000);
    std::cout << mpfr_float::default_precision() << std::endl;
    std::cout << sqrt(a) << std::endl; // print root-2

    // Operations at fixed precision and full numeric_limits support:
    mpfr_float_100 b = 2;
    std::cout << std::numeric_limits<mpfr_float_100>::digits << std::endl;
    // We can use any C++ std lib function:
    std::cout << log(b) << std::endl; // print log(2)
    // We can also use any function from Boost.Math:
    std::cout << boost::math::tgamma(b) << std::endl;
    // These even work when the argument is an expression template:
    std::cout << boost::math::tgamma(b * b) << std::endl;

    // Access the underlying data:
    mpfr_t r;
    mpfr_init(r);
    mpfr_set(r, b.backend().data(), GMP_RNDN);
    mpfr_clear(r);

    //CppAD::AD<mpfr_float> ca = mpfr_float(0.0);
    //CppAD::AD<mpfr_float> cb = mpfr_float(1.0);
    //CppAD::AD<mpfr_float> dd = ca + cb;
    mpfr_float ta = mpfr_float(1.0);
    cout << ta.is_zero() << endl;

    vector<CppAD::AD<mpfr_float>> s_k(1);
    vector<mpfr_float> s_k_x(1);
    s_k[0] = mpfr_float(1.0);
    s_k_x[0] = mpfr_float(1.0);
    CppAD::Independent(s_k);

    ADmpfr_float aa = mpfr_float(5.0);
    ADmpfr_float bb = mpfr_float(20.0);
    vector<ADmpfr_float> cc(1);
    cc[0] = aa * s_k[0];
    cc[0] = cc[0] * bb * bb;
    cout << cc[0] << endl;
    CppAD::ADFun<mpfr_float> N_s(s_k, cc);
    vector<mpfr_float> jac_N(1);
    jac_N = N_s.Jacobian(s_k_x);
    cout << jac_N[0] << endl;
    return 0;
}
//int main() {
    //mp::mpfr_float::default_precision(10000);
    //std::cout.precision(100);
    ////mp::mpfr_float_1000 c = 0.123456789123456789123456789;
    //mp::mpfr_float c("0.123456789123456789123456789");
    //std::cout << "Value of c: " << c << std::endl;
    //// Define two high-precision numbers
    //mp::mpfr_float a("123456789.123456789123456789");
    //mp::mpfr_float b("987654321.987654321987654321");
    //// Perform arithmetic operations
    //mp::mpfr_float sum = a + b;
    //mp::mpfr_float product = a * b;
    //// Set the precision for output to 50 digits
    //// Output the results
    //std::cout << "Sum of a and b: " << sum << std::endl;
    //std::cout << "Product of a and b: " << product << std::endl;
    //double a1 = 123456789.123456789123456789;
    //double b1 = 987654321.987654321987654321;
    //// Perform arithmetic operations
    //double sum1 = a1 + b1;
    //double product1 = a1 * b1;
    //// Set the precision for output to 50 digits
    //std::cout.precision(50);
    //// Output the results
    //std::cout << "Sum of a and b: " << sum1 << std::endl;
    //std::cout << "Product of a and b: " << product1 << std::endl;
    //return 0;
//}


//#include <stdio.h>
//#include <gmp.h>
//
//int main() {
//    mpz_t a, b, sum, difference, product, quotient;
//    mpz_init(a);
//    mpz_init(b);
//    mpz_init(sum);
//    mpz_init(difference);
//    mpz_init(product);
//    mpz_init(quotient);
//
//    mpz_set_str(a, "123456789012345678901234567890", 10); // 设置a为一个大整数
//    mpz_set_str(b, "987654321098765432109876543210", 10); // 设置b为一个大整数
//
//    // 计算和
//    mpz_add(sum, a, b);
//
//    // 计算差
//    mpz_sub(difference, a, b);
//
//    // 计算乘积
//    mpz_mul(product, a, b);
//
//    // 计算除法
//    mpz_div(quotient, a, b);
//
//    gmp_printf("a + b = %Zd\n", sum);
//    gmp_printf("a - b = %Zd\n", difference);
//    gmp_printf("a * b = %Zd\n", product);
//    gmp_printf("a / b = %Zd\n", quotient);
//
//    mpz_clear(a);
//    mpz_clear(b);
//    mpz_clear(sum);
//    mpz_clear(difference);
//    mpz_clear(product);
//    mpz_clear(quotient);
//
//    return 0;
//}
