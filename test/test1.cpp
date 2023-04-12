#include "stdio.h"
#include "stdlib.h"
#include "iostream"
#include "../utils/util.hpp"
#include "../utils/PRgenerator.hpp"
#include "omp.h"
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
void test1()
{
    srand(time(NULL));
    for (size_t i = 0; i < 10; i++)
    {
        size_t cols = rand() % 10;
        size_t rows = rand() % 10;
        MatrixXl rand_mtx = rand_matrixXl(rows, cols);
        string temp;
        matrix2string(rand_mtx, temp);
        MatrixXl _rand_mtx(rows, cols);
        string2matrix(temp, _rand_mtx);
        assert(rand_mtx == _rand_mtx);
    }
    for (size_t i = 0; i < 10; i++)
    {
        size_t vec_size = rand() % 10000;
        vector<int64_t> vec = rand_vec_int64(vec_size);
        string temp;
        vector2string(vec, temp);
        vector<int64_t> _vec;
        string2vector(temp, _vec);
        assert(vec == _vec);
    }
}
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
    }// implict barrier

    printf("all thread have done its work\n");
}

int main(int argc, char const *argv[])
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
    test1();
    test2();
    test3();
}