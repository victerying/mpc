#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "../utils/PRgenerator.hpp"
#include "../utils/SocketMessenger.hpp"
#include "../utils/util.hpp"
#include "../mpc_engine/ABY3Engine.hpp"

/// @brief 测试share和reveal (vecI64, matrix, vecbools)
/// @param aby3Engine
void test1(ABY3Engine &aby3Engine)
{
    size_t nextParty = (aby3Engine.party_id + 1) % 3;
    size_t prevParty = (aby3Engine.party_id + 2) % 3;
    PRgenerator commRando(20, 0);
    // VecI64
    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64();
        }
        array<VecI64, 2> Xj;
        aby3Engine.share(X_plain, Xj);
        // printf("X_plain\n");
        // printVecI64(X_plain);
        // printf("Xj[0]\n");
        // printVecI64(Xj[0]);
        // printf("Xj[1]\n");
        // printVecI64(Xj[1]);
        VecI64 _X_plain;
        aby3Engine.reveal(Xj, _X_plain);
        // printf("_X_plain\n");
        // printVecI64(_X_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(X_plain[i] == _X_plain[i]);
        }
    }
    // MatrixXl
    {
        size_t rows = 5, cols = 4;
        MatrixXl X_plain(rows, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                X_plain(i, j) = commRando.pop_int64();
            }
        }

        array<MatrixXl, 2> Xj;

        aby3Engine.share_matrix(X_plain, Xj);
        // printf("X_plain\n");
        // std::cout << X_plain << "\n";
        // printf("Xj[0]\n");
        // std::cout << Xj[0] << "\n";
        // printf("Xj[1]\n");
        // std::cout << Xj[1] << "\n";

        MatrixXl _X_plain(rows, cols);
        aby3Engine.reveal_matrix(Xj, _X_plain);
        // printf("_X_plain\n");
        // std::cout << _X_plain << "\n";
        assert(X_plain == _X_plain);
    }
    // vecbool
    {
        size_t vec_size = 10;
        vecbool X_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_bool();
        }
        array<vecbool, 2> Xj;
        aby3Engine.share_bools(X_plain, Xj);
        vecbool _X_plain;
        aby3Engine.reveal_bools(Xj, _X_plain);
        // printf("X_plain\n");
        // printVecbool(X_plain);
        // printf("Xj[0]\n");
        // printVecbool(Xj[0]);
        // printf("Xj[1]\n");
        // printVecbool(Xj[1]);
        // printf("_X_plain\n");
        // printVecbool(_X_plain);

        assert(X_plain == _X_plain);
    }
}

/// @brief 测试乘法
/// @param aby3Engine
void test2(ABY3Engine &aby3Engine)
{
    size_t nextParty = (aby3Engine.party_id + 1) % 3;
    size_t prevParty = (aby3Engine.party_id + 2) % 3;
    PRgenerator commRando(time(NULL), 0);

    {
        size_t rows = 5, lens = 4, cols = 3;
        MatrixXl X_plain(rows, lens), Y_plain(lens, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < lens; j++)
            {
                X_plain(i, j) = commRando.pop_int64();
            }
        }
        for (size_t i = 0; i < lens; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                Y_plain(i, j) = commRando.pop_int64();
            }
        }
        array<MatrixXl, 2> Xj, Yj, Zj;
        aby3Engine.share_matrix(X_plain, Xj);
        aby3Engine.share_matrix(Y_plain, Yj);
        aby3Engine.matrix_mul(Xj, Yj, Zj);
        MatrixXl Z_plain;
        aby3Engine.reveal_matrix(Zj, Z_plain);
        // printf("Z_plain\n");
        // std::cout << Z_plain << "\n";
        // printf("X_plain * Y_plain\n");
        // std::cout << X_plain * Y_plain << "\n";
        assert(Z_plain == X_plain * Y_plain);
    }

    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size), Y_plain(vec_size), Z_plain;
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64();
            Y_plain[i] = commRando.pop_int64();
        }
        array<VecI64, 2> Xj, Yj, Zj;
        aby3Engine.share(X_plain, Xj);
        aby3Engine.share(Y_plain, Yj);
        aby3Engine.vector_mul(Xj, Yj, Zj);
        aby3Engine.reveal(Zj, Z_plain);
        VecI64 temp = X_plain * Y_plain;
        // printf("Z_plain\n");
        // printVecI64(Z_plain);
        // printf("X_plain *Y_plain \n");
        // printVecI64(temp );
        for (size_t i = 0; i < vec_size; i++)
        {
            temp[i] = Z_plain[i];
        }
    }

    {
        size_t vec_size = 100;
        vecbool X_plain(vec_size), Y_plain(vec_size), Z_plain;
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_bool();
            Y_plain[i] = commRando.pop_bool();
        }
        array<vecbool, 2> Xj, Yj, Zj;
        aby3Engine.share_bools(X_plain, Xj);
        aby3Engine.share_bools(Y_plain, Yj);
        aby3Engine.vecbool_mul(Xj, Yj, Zj);
        aby3Engine.reveal_bools(Zj, Z_plain);
        // printf("X_plain\n");
        // printVecbool(X_plain);
        // printf("Y_plain\n");
        // printVecbool(Y_plain);
        // printf("Z_plain\n");
        // printVecbool(Z_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(X_plain[i] * Y_plain[i] == Z_plain[i]);
        }
    }
}

/// @brief 测试msb
/// @param aby3Engine
void test3(ABY3Engine &aby3Engine)
{
    size_t nextParty = (aby3Engine.party_id + 1) % 3;
    size_t prevParty = (aby3Engine.party_id + 2) % 3;
    PRgenerator commRando(time(NULL), 0);
    {
        size_t vec_size = 4;
        VecI64 X_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64();
        }
        array<VecI64, 2> Xj, Gj, Pj;
        aby3Engine.share(X_plain, Xj);
        aby3Engine.add_prepare(Xj, Gj, Pj);
        array<vecbool, 2> Zj;
        aby3Engine.msb(Gj, Pj, Zj);
        vecbool Z_debug;
        aby3Engine.reveal_bools(Zj, Z_debug);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(Z_debug[i] == (X_plain[i] < 0));
        }
        array<VecI64, 2> Wj;
        VecI64 W_debug;
        aby3Engine.bit_inject(Xj, Zj, Wj);
        aby3Engine.reveal(Wj, W_debug);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(X_plain[i] * Z_debug[i] == W_debug[i]);
        }
        
        // printf("X_plain:\n");
        // printVecI64(X_plain);
        // printf("Z_debug:\n");
        // printVecbool(Z_debug);
        // printf("W_debug:\n");
        // printVecI64(W_debug);
    }
}

/// @brief 测试trunc
/// @param aby3Engine
void test4(ABY3Engine &aby3Engine)
{
    size_t nextParty = (aby3Engine.party_id + 1) % 3;
    size_t prevParty = (aby3Engine.party_id + 2) % 3;
    PRgenerator commRando(time(NULL), 0);
    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size), Y_plain;
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64();
            X_plain[i] = X_plain[i] >> 32;
        }
        array<VecI64, 2> Xj, Yj;
        aby3Engine.share(X_plain, Xj);
        aby3Engine.trunc(Xj, Yj, 16);
        aby3Engine.reveal(Yj, Y_plain);
        // printf("X_plain:\n");
        // printVecI64(X_plain);
        // printf("Y_plain:\n");
        // printVecI64(Y_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            int64_t temp = X_plain[i] >> 16;
            assert(Y_plain[i] == temp || Y_plain[i] == (temp - 1));
        }
    }
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
    ABY3Engine aby3Engine(&socketMessenger, party_id);
    test1(aby3Engine);
    test2(aby3Engine);
    test3(aby3Engine);
    test4(aby3Engine);
    return 0;
}
