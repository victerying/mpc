#pragma once
#include "../eigen-3.4.0/Eigen/Dense"
#include <vector>
using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
using std::vector;
/// @brief all of parameters are secret share
class MpcEngine
{
private:
public:
    MpcEngine() = default;

    /// @brief c = a * b (matrix multiply over ring Z 2^64),
    /// @param a input MatrixXl
    /// @param b input MatrixXl
    /// @param c output MatrixXl
    virtual void matrix_mul(const MatrixXl &a, const MatrixXl &b, MatrixXl &c) = 0;

    /// @brief c[i] = a[i] * b[i] (element wise multiply over ring Z 2^64)
    /// @param a
    /// @param b
    /// @param c
    virtual void dot_mul(const vector<int64_t> &a, const vector<int64_t> &b, vector<int64_t> &c) = 0;

    /// @brief a[i] = a[i] >> precision (element wise shft left)
    /// @param a input output
    /// @param precision public known
    virtual void dot_shft_left(vector<int64_t> &a, size_t precision) = 0;

    /// @brief a(i,j) = a(i,j)  >> precision ( element wise shft left)
    /// @param a input output
    /// @param precision  public known
    virtual void matrix_shft_left(MatrixXl &a, size_t precision) = 0;

    /// @brief b(i,j) = msb{ a(i,j) } msb stands for most significant bit
    /// @param a input
    /// @param b output
    virtual void matrix_msb(const MatrixXl &a, MatrixXl &b) = 0;


    /// @brief b[i]= msb{ a[i] } msb stands for most significant bit
    /// @param a input
    /// @param b output
    virtual void dot_msb(const vector<int64_t> &a, vector<int64_t> &b) = 0;


    ~MpcEngine() = default;
};