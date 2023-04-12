#include "../eigen-3.4.0/Eigen/Dense"
#include <stdlib.h>
#include <string>
#include <vector>
using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
using std::string;
using std::vector;
static void matrix2string(const MatrixXl &input, string &output)
{
    size_t rows = input.rows();
    size_t cols = input.cols();
    output.resize(rows * cols * sizeof(int64_t));
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            memcpy((char *)output.c_str() + (i * cols + j) * sizeof(int64_t), &input(i, j), sizeof(int64_t));
        }
    }
}

/// @brief before call ,the rows and cols of output must be set
/// @param input
/// @param output
static void string2matrix(const string &input, MatrixXl &output)
{
    size_t rows = output.rows();
    size_t cols = output.cols();
    assert(input.size() == rows * cols * sizeof(int64_t));
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            output(i, j) = *(int64_t *)(input.c_str() + (i * cols + j) * sizeof(int64_t));
        }
    }
}

static void vector2string(const vector<int64_t> &input, string &output)
{
    size_t vec_size = input.size();
    output.resize(vec_size * sizeof(int64_t));
    for (size_t i = 0; i < vec_size; i++)
    {
        memcpy((char *)output.c_str() + i * sizeof(int64_t), &input[i], sizeof(int64_t));
    }
    
}
static void string2vector(const string &input, vector<int64_t> &output)
{
    assert(input.size() %sizeof(int64_t) == 0);
    size_t vec_size = input.size() / sizeof(int64_t);
    output.resize(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        output[i] = *(int64_t *)(input.c_str() + i * sizeof(int64_t) );
    }
}