#include "stdio.h"
#include "stdlib.h"
#include "iostream"
#include "../utils/util.hpp"
#include "../utils/PRgenerator.hpp"
#include "omp.h"
using zpshare = std::array<uint16_t, 64>;
static inline int64_t rand_int64()
{
    return (int64_t)rand() << 62 | (int64_t)rand() << 31 | (int64_t)rand();
}
static inline vector<int64_t> rand_vec_int64(size_t vec_size)
{
    vector<int64_t> ret(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        ret[i] = rand_int64();
    }
    return ret;
}
static inline MatrixXl rand_matrixXl(size_t rows, size_t cols)
{
    MatrixXl ret(rows, cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            ret(i, j) = rand_int64();
        }
    }
    return ret;
}
static inline vector<zpshare> rand_zpshare(size_t vec_size)
{
    vector<zpshare> ret(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        for (size_t j = 0; j < 64; j++)
        {
            ret[i][j] = rand();
        }
    }
    return ret;
}

static const uint64_t L_1 = -1;

static inline uint64_t add_L_1(uint64_t a, uint64_t b)
{
    bool carry;
    uint64_t c = a + b;
    carry = c < a;
    return (c + carry) % L_1;
}

static inline uint64_t sub_L_1(uint64_t a, uint64_t b)
{

    uint64_t _b = (L_1 - b) % L_1;
    return add_L_1(a, _b);
}

/// @brief 测试c语言基础
void test0()
{

    assert(sizeof(size_t) == 8);
    assert(sizeof(int) == 4);
    assert(sizeof(short) == 2);
    assert(sizeof(char) == 1);

    assert(INT32_MAX == 0x7fffffff);
    char buf[100];
    size_t shift = 50;
    void *ptr1 = buf;
    void *ptr2 = ptr1 + shift;
    assert(ptr1 == &buf[0]);
    assert(ptr2 == &buf[shift]);

    string str1, str2;
    str1.resize(rand() % 100);
    str2.resize(rand() % 100);
    assert((str1 + str2).size() == str1.size() + str2.size());

    for (size_t i = 0; i < 10; i++)
    {
        int64_t a = rand_int64();
        int64_t b = rand_int64();
        int64_t c = rand_int64();
        int64_t r1 = (a + b) * c;
        int64_t r2 = a * c + b * c;
        assert(r1 == r2);
    }
}

/// @brief 测试util中各种序列化
void test1()
{
    srand(time(NULL));
    for (size_t i = 0; i < 10; i++)
    {
        size_t cols = rand() % 10 + 5;
        size_t rows = rand() % 10 + 5;
        MatrixXl rand_mtx = rand_matrixXl(rows, cols);
        string temp;
        matrix2string(rand_mtx, temp);
        MatrixXl _rand_mtx(rows, cols);
        string2matrix(temp, _rand_mtx);
        assert(rand_mtx == _rand_mtx);
    }
    for (size_t i = 0; i < 10; i++)
    {
        size_t vec_size = rand() % 1000 + 500;
        vector<int64_t> vec = rand_vec_int64(vec_size);
        string temp;
        vector2string(vec, temp);
        vector<int64_t> _vec;
        string2vector(temp, _vec);
        assert(vec == _vec);
    }
    for (size_t i = 0; i < 10; i++)
    {
        size_t vec_size = rand() % 100 + 100;
        vector<zpshare> vec_zpshare = rand_zpshare(vec_size);
        string temp;
        vecZpshare2string(vec_zpshare, temp);
        vector<zpshare> _vec_zpshare;
        string2vecZpshare(temp, _vec_zpshare);
        assert(vec_zpshare == _vec_zpshare);
        if (i == 0)
        {
            for (size_t k = 0; k < 12; k++)
            {
                printf("%d ", vec_zpshare[50][k + 30]);
            }
            printf("\n");
            for (size_t k = 0; k < 12; k++)
            {
                printf("%d ", _vec_zpshare[50][k + 30]);
            }
            printf("\n");
        }
    }
    for (size_t i = 0; i < 10; i++)
    {
        size_t cols = rand() % 10 + 5;
        size_t rows = rand() % 10 + 5;
        MatrixXl rand_mtx = rand_matrixXl(rows, cols);
        vector<int64_t> temp;
        matrix2vector(rand_mtx, temp);
        MatrixXl _rand_mtx(rows, cols);
        vector2matrix(temp, _rand_mtx);
        assert(rand_mtx == _rand_mtx);
    }
}

/// @brief 测试伪随机数生成器
void test2()
{
    uint64_t seed = rand_int64();
    uint64_t lcg = rand_int64();
    PRgenerator pr_generator(seed, lcg);
    // for (size_t i = 0; i < 10; i++)
    // {
    //     for (size_t j = 0; j < 5; j++)
    //     {
    //         printf("%lx ", pr_generator.int_buff[i * 5 + j]);
    //     }
    //     printf("\n");
    // }
    // for (size_t i = 0; i < 10; i++)
    // {
    //     for (size_t j = 0; j < 5; j++)
    //     {
    //         printf("%d ", pr_generator.bool_buff[i * 5 + j]);
    //     }
    //     printf("\n");
    // }
}

/// @brief 测试openmp
void test3()
{
#pragma omp parallel num_threads(2)
#pragma omp sections
    {
#pragma omp section
        {
            int thread_id = omp_get_thread_num();
            int num_threads = omp_get_num_threads();
            printf("there are %d threads, thread %d is doing section 1\n", num_threads, thread_id);
        }
#pragma omp section
        {
            int thread_id = omp_get_thread_num();
            int num_threads = omp_get_num_threads();
            printf("there are %d threads, thread %d is doing section 2\n", num_threads, thread_id);
        }
    } // implict barrier

    printf("all thread have done its work\n");
}

/// @brief 测试L-1上的加法减法
void test4()
{
    printf("test4\n");
    srand(time(NULL));
    for (size_t i = 0; i < 100; i++)
    {
        {
            int8_t flag = rand() % 2;
            uint64_t temp = flag;
            int64_t add1 = rand_int64();
            int64_t add2 = sub_L_1(temp, add1);
            uint64_t _temp = add_L_1(add1, add2);
            assert(temp == _temp);
        }
        {
            uint64_t temp = rand_int64();
            int64_t add1 = rand_int64();
            int64_t add2 = sub_L_1(temp, add1);
            uint64_t _temp = add_L_1(add1, add2);
            assert(temp == _temp);
        }
    }
}

int main(int argc, char const *argv[])
{
    printf("%lu\n", sizeof(long long));
    test0();
    test1();
    test2();
    test3();
    test4();
}