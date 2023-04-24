#include <iostream>
#include "stdlib.h"
#include "stdio.h"
#include "../utils/PRgenerator.hpp"
#include "../utils/SocketMessenger.hpp"
#include "../utils/util.hpp"
#include "../mpc_engine/ABY3Engine.hpp"
#include "../mpc_engine/MyEngine.hpp"
#include "../mpc_engine/SecurennEngine.hpp"

void bench_aby3(SocketMessenger &socketMessenger)
{
    ABY3Engine aby3Engine(&socketMessenger, socketMessenger.party_id);
    size_t repeat_time = 10;
    PRgenerator commRando(time(NULL), 0);
    size_t width = 200;
    clock_t spend = 0;
    for (size_t t = 0; t < repeat_time; t++)
    {
        MatrixXl X_plain(width, width), Y_plain(width, width);
        for (size_t i = 0; i < width; i++)
        {
            for (size_t j = 0; j < width; j++)
            {
                X_plain(i, j) = commRando.pop_int64();
                Y_plain(i, j) = commRando.pop_int64();
            }
        }
        array<MatrixXl, 2> Xj, Yj, Zj;
        aby3Engine.share_matrix(X_plain, Xj);
        aby3Engine.share_matrix(Y_plain, Yj);

        clock_t start1 = clock();
        aby3Engine.matrix_mul(Xj, Yj, Zj);
        clock_t end1 = clock();

        Xj = std::move(Zj);
        array<VecI64, 2> Xj_vec, Zj_vec;
        for (size_t i = 0; i < 2; i++)
            matrix2vecI64(Xj[i], Xj_vec[i]);

        clock_t start2 = clock();
        aby3Engine.trunc(Xj_vec, Zj_vec, 16);
        clock_t end2 = clock();

        for (size_t i = 0; i < 2; i++)
            vecI642matrix(Zj_vec[i], Zj[i]);

        spend += (end1 - start1);
        spend += (end2 - start2);
    }
    spend = spend / repeat_time;
    printf("[LOG]:\tbench_aby3 matrix_mul & trunction\n");
    print_plus();
    printf("spent time: %ld\n", spend);
    print_minus();

    spend = 0;
    for (size_t t = 0; t < repeat_time; t++)
    {
        size_t vec_size = width * 10;
        VecI64 X_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64();
        }
        array<VecI64, 2> Xj, Gj, Pj;
        array<vecbool, 2> Zj;
        aby3Engine.share(X_plain, Xj);
        clock_t start = clock();
        aby3Engine.add_prepare(Xj, Gj, Pj);
        aby3Engine.msb(Gj, Pj, Zj);
        aby3Engine.bit_inject(Xj, Zj, Gj);
        clock_t end = clock();
        spend += (end - start);
    }
    spend = spend / repeat_time;
    printf("[LOG]:\tbench_aby3 msb & bit_inject\n");
    print_plus();
    printf("spent time: %ld\n", spend);
    print_minus();
}

void bench_securenn(SocketMessenger &socketMessenger)
{
    SecurennEngine securenn_eigine(&socketMessenger, socketMessenger.party_id);
    size_t repeat_time = 10;
    PRgenerator commRando(time(NULL), 0);
    size_t width = 200;
    clock_t spend = 0;
    for (size_t t = 0; t < repeat_time; t++)
    {
        MatrixXl Xj(width, width), Yj(width, width), Zj(width, width);
        vector<int64_t> X_plain_vec(width * width), Y_plain_vec(width * width), Xj_vec(width * width), Yj_vec(width * width);
        for (size_t i = 0; i < width * width; i++)
        {
            X_plain_vec[i] = commRando.pop_int64() / 2;
            Y_plain_vec[i] = commRando.pop_int64() / 2;
        }
        securenn_eigine.share(X_plain_vec, Xj_vec);
        securenn_eigine.share(Y_plain_vec, Yj_vec);
        vector2matrix(Xj_vec, Xj);
        vector2matrix(Yj_vec, Yj);
        clock_t start = clock();
        securenn_eigine.matrix_mul(Xj, Yj, Zj);
        for (size_t i = 0; i < width; i++)
        {
            for (size_t j = 0; j < width; j++)
            {
                Zj(i, j) = Zj(i, j) >> 16;
            }
        }
        clock_t end = clock();
        spend += (end - start);
    }
    spend = spend / repeat_time;
    printf("[LOG]:\tbench_securenn matrix_mul & trunction\n");
    print_plus();
    printf("spent time: %ld\n", spend);
    print_minus();

    spend = 0;
    for (size_t t = 0; t < repeat_time; t++)
    {
        size_t vec_size = width * 10;
        vector<int64_t> A_plain(vec_size), Aj, Bj, Cj;
        for (size_t i = 0; i < vec_size; i++)
        {
            A_plain[i] = commRando.pop_int64();
            A_plain[i] = A_plain[i] / 2;
        }
        securenn_eigine.share(A_plain, Aj);

        clock_t start = clock();
        securenn_eigine.msb(Aj, Bj);
        securenn_eigine.vector_mul(Aj, Bj, Cj);
        clock_t end = clock();
        spend += (end - start);
    }
    spend = spend / repeat_time;
    printf("[LOG]:\tbench_securenn msb & bit_inject\n");
    print_plus();
    printf("spent time: %ld\n", spend);
    print_minus();
}

void bench_me(SocketMessenger &socketMessenger)
{
    MyEngine myengine(&socketMessenger, socketMessenger.party_id);
    size_t repeat_time = 10;
    PRgenerator commRando(time(NULL), 0);
    size_t width = 200;
    clock_t spend = 0;
    for (size_t t = 0; t < repeat_time; t++)
    {
        MatrixXl X_plain(width, width), Y_plain(width, width), Z_plain(width, width);
        MatrixXl Xj, Yj, Zj;
        for (size_t i = 0; i < width; i++)
            for (size_t k = 0; k < width; k++)
            {
                X_plain(i, k) = commRando.pop_int64();
                Y_plain(i, k) = commRando.pop_int64();
            }
        myengine.share_matrix(X_plain, Xj);
        myengine.share_matrix(Y_plain, Yj);
        clock_t start = clock();
        myengine.matrix_mul(Xj, Yj, Zj);
        for (size_t i = 0; i < width; i++)
        {
            for (size_t j = 0; j < width; j++)
            {
                Zj(i, j) = Zj(i, j) >> 16;
            }
        }
        clock_t end = clock();
        spend += (end - start);
    }
    spend = spend / repeat_time;
    printf("[LOG]:\tbench_me matrix_mul & trunction\n");
    print_plus();
    printf("spent time: %ld\n", spend);
    print_minus();

    spend = 0;
    for (size_t t = 0; t < repeat_time; t++)
    {
        size_t vec_size = width * 10;
        VecI64 X_plain(vec_size), Xj(vec_size);
        vecbool Aj(vec_size);
        VecI64 Zj(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64() / 2;
        }
        myengine.share(X_plain, Xj);

        clock_t start = clock();
        myengine.msb(Xj, Aj);
        myengine.bit_inject(Xj, Aj, Zj);
        clock_t end = clock();
        spend += (end - start);
    }
    spend = spend / repeat_time;
    printf("[LOG]:\tbench_me msb & bit_inject\n");
    print_plus();
    printf("spent time: %ld\n", spend);
    print_minus();
}

int main(int argc, char const *argv[])
{
    if (argc != 2)
    {
        printf("usage:\ntest2 <partyid>\n");
        return 1;
    }
    size_t party_id = strtol(argv[1], NULL, 10);
    if (party_id > 2)
    {
        printf("party_id must be in {1,2,3}\n");
        return 1;
    }
    SocketMessenger socketMessenger(party_id);
    printf("[INFO]:\tCLOCKS_PER_SEC:%ld\n", CLOCKS_PER_SEC);

    bench_aby3(socketMessenger);
    bench_securenn(socketMessenger);
    bench_me(socketMessenger);
    return 0;
}
