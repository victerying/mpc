#pragma once
#include "../eigen-3.4.0/Eigen/Dense"
#include <stdlib.h>
#include <string>
#include <vector>
#include <array>
using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
using VecI64 = Eigen::Array<int64_t, Eigen::Dynamic, 1>;
using zpshare = std::array<uint16_t, 64>;
using vecbool = vector<int8_t>;

using std::string;
using std::vector;

/// @brief convert matrix to string
/// @param input (rows * cols)
/// @param output size = rows * cols * sizeof(int64_t) (after call)
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

/// @brief convert string to matrix
/// @param input size = rows * cols * sizeof(int64_t)
/// @param output before call, the rows and cols of output must be set correctly
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

/// @brief convert vector<int64_t> to string
/// @param input
/// @param output str_size = vec_size * sizeof(int64_t) (after call)
static void vector2string(const vector<int64_t> &input, string &output)
{
    size_t vec_size = input.size();
    output.resize(vec_size * sizeof(int64_t));
    for (size_t i = 0; i < vec_size; i++)
    {
        memcpy((char *)output.c_str() + i * sizeof(int64_t), &input[i], sizeof(int64_t));
    }
}

/// @brief convert string to vector<int64_t>
/// @param input str_size %8=0  (8 = sizeof(int64_t))
/// @param output vec_size = str_size /8
static void string2vector(const string &input, vector<int64_t> &output)
{
    assert(input.size() % sizeof(int64_t) == 0);
    size_t vec_size = input.size() / sizeof(int64_t);
    output.resize(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        output[i] = *(int64_t *)(input.c_str() + i * sizeof(int64_t));
    }
}

static void vecZpshare2string(const vector<zpshare> &input, string &output)
{
    size_t vec_size = input.size();
    size_t element_size = 64 * sizeof(uint16_t);
    output.resize(vec_size * element_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        for (size_t j = 0; j < 64; j++)
        {
            memcpy((char *)output.c_str() + i * 64 * 2 + j * 2, &input[i][j], 2);
        }
    }
}

static void string2vecZpshare(const string &input, vector<zpshare> &output)
{
    size_t element_size = 64 * sizeof(uint16_t);
    size_t str_size = input.size();
    assert(str_size % element_size == 0);
    size_t vec_size = str_size / element_size;
    output.resize(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        for (size_t j = 0; j < 64; j++)
        {
            memcpy(&output[i][j], input.c_str() + i * 64 * 2 + j * 2, 2);
        }
    }
}

static void matrix2vector(const MatrixXl &input, vector<int64_t> &output)
{
    size_t rows = input.rows();
    size_t cols = input.cols();
    output.resize(rows * cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            output[i * cols + j] = input(i, j);
        }
    }
}

static void vector2matrix(const vector<int64_t> &input, MatrixXl &output)
{
    size_t rows = output.rows();
    size_t cols = output.cols();
    assert(input.size() == rows * cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            output(i, j) = input[i * cols + j];
        }
    }
}

static void vecI642string(const VecI64 &input, string &output)
{
    size_t vec_size = input.size();
    output.resize(vec_size * sizeof(int64_t));
    for (size_t i = 0; i < vec_size; i++)
    {
        memcpy((char *)output.c_str() + i * sizeof(int64_t), &input(i), sizeof(int64_t));
    }
}

static void string2vecI64(const string &input, VecI64 &output)
{
    assert(input.size() % sizeof(int64_t) == 0);
    size_t vec_size = input.size() / sizeof(int64_t);
    output.resize(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        output(i) = *(int64_t *)(input.c_str() + i * sizeof(int64_t));
    }
}

static void vecI642matrix(const VecI64 &input, MatrixXl &output)
{
    size_t rows = output.rows();
    size_t cols = output.cols();
    assert(input.size() == rows * cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t k = 0; k < cols; k++)
        {
            output(i, k) = input[i * cols + k];
        }
    }
}

static void matrix2vecI64(const MatrixXl &input, VecI64 &output)
{
    size_t rows = input.rows();
    size_t cols = input.cols();
    output.resize(rows * cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t k = 0; k < cols; k++)
        {
            output[i * cols + k] = input(i, k);
        }
    }
}

static void vecbool2string(const vector<int8_t> &input, string &output)
{
    output.resize(input.size());
    for (size_t i = 0; i < input.size(); i++)
    {
        output[i] = input[i];
    }
}

static void string2vecbool(const string &input, vector<int8_t> &output)
{
    output.resize(input.size());
    for (size_t i = 0; i < input.size(); i++)
    {
        output[i] = input[i];
    }
}

static void print_vec_int64(const vector<int64_t> &input)
{
    for (size_t i = 0; i < input.size(); i++)
    {
        printf("%016lX ", input[i]);
    }
    printf("\n");
}

static void print_vec_zpshare(const vector<zpshare> &input)
{
    for (size_t i = 0; i < input.size(); i++)
    {
        for (size_t k = 0; k < 64; k++)
        {
            printf("%2d ", input[i][k]);
            if (k % 32 == 31)
                printf("\n");
        }
        printf("\n");
    }
}

static void printVecI64(const VecI64 &input)
{
    for (size_t i = 0; i < input.size(); i++)
    {
        printf("%016lX ", input[i]);
    }
    printf("\n");
}

static void printVecbool(const vecbool &input)
{
    for (size_t i = 0; i < input.size(); i++)
    {
        printf("%d ", input[i]);
    }
    printf("\n");
}
