#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "../utils/PRgenerator.hpp"
#include "../utils/SocketMessenger.hpp"
#include "../utils/util.hpp"
#include "../mpc_engine/MyEngine.hpp"

void test1(MyEngine &myengine)
{

    PRgenerator commRando(time(NULL), 0);
    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
            X_plain[i] = commRando.pop_int64();

        VecI64 Xj, _X_plain;
        myengine.share(X_plain, Xj);
        myengine.reveal(Xj, _X_plain);
        for (size_t i = 0; i < vec_size; i++)
            assert(X_plain[i] == _X_plain[i]);

        // printf("X_plain\n");
        // printVecI64(X_plain);
        // printf("_X_plain\n");
        // printVecI64(_X_plain);
    }
    {
        size_t rows = 9, cols = 10;
        MatrixXl X_plain(rows, cols);
        MatrixXl Xj, _X_plain;
        for (size_t i = 0; i < rows; i++)
            for (size_t k = 0; k < cols; k++)
                X_plain(i, k) = commRando.pop_int64();
        myengine.share_matrix(X_plain, Xj);
        myengine.reveal_matrix(Xj, _X_plain);
        assert(X_plain == _X_plain);
        // printf("X_plain\n");
        // std::cout << X_plain << "\n";
        // printf("_X_plain\n");
        // std::cout << _X_plain << "\n";
    }
    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size), Y_plain(vec_size), Z_plain, Xj, Yj, Zj;
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64();
            Y_plain[i] = commRando.pop_int64();
        }
        myengine.share(X_plain, Xj);
        myengine.share(Y_plain, Yj);
        myengine.vector_mul(Xj, Yj, Zj);
        myengine.reveal(Zj, Z_plain);
        VecI64 temp_mul = X_plain * Y_plain;
        for (size_t i = 0; i < vec_size; i++)
            assert(Z_plain[i] == temp_mul[i]);

        // printf("Z_plain\n");
        // printVecI64(Z_plain);
        // printf("X*Y\n");
        // printVecI64(X_plain * Y_plain);
    }
    {
        size_t rows = 9, lens = 11, cols = 10;
        MatrixXl X_plain(rows, lens), Y_plain(lens, cols), Z_plain(rows, cols);
        MatrixXl Xj, Yj, Zj;

        for (size_t i = 0; i < rows; i++)
            for (size_t k = 0; k < lens; k++)
                X_plain(i, k) = commRando.pop_int64();
        for (size_t i = 0; i < lens; i++)
            for (size_t k = 0; k < cols; k++)
                Y_plain(i, k) = commRando.pop_int64();
        myengine.share_matrix(X_plain, Xj);
        myengine.share_matrix(Y_plain, Yj);
        myengine.matrix_mul(Xj, Yj, Zj);
        myengine.reveal_matrix(Zj, Z_plain);
        assert(Z_plain == X_plain * Y_plain);
        // printf("Z_plain\n");
        // std::cout << Z_plain << "\n";
        // printf("X*Y\n");
        // std::cout << X_plain * Y_plain << "\n";
    }
}

void test2(MyEngine &myengine)
{
    PRgenerator commRando(time(NULL), 0);

    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size), R_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64() % L_1;
            R_plain[i] = commRando.pop_int64() % L_1;
        }
        {
            X_plain[0] = 0;
            X_plain[1] = -2;
            R_plain[3] = 0;
            R_plain[4] = -2;
            X_plain[5] = R_plain[5];
            X_plain[6] = (R_plain[6] + 1) % L_1;
            X_plain[7] = (R_plain[7] - 1) % L_1;
        }
        vector<zpshare> Xbits_plain(vec_size), X_bits[2];
        X_bits[0].resize(vec_size);
        X_bits[1].resize(vec_size);
        for (size_t k = 0; k < vec_size; k++)
        {
            bitset<64> btst_x(X_plain[k]);
            for (size_t l = 0; l < 64; l++)
            {
                Xbits_plain[k][l] = (bool)btst_x[l];
                X_bits[0][k][l] = commRando.pop_int64();
                X_bits[0][k][l] = X_bits[0][k][l] % P;

                X_bits[1][k][l] = Xbits_plain[k][l] + P - X_bits[0][k][l];
                X_bits[1][k][l] = X_bits[1][k][l] % P;
            }
        }
        vecbool Beta1_plain(vec_size), Beta2_plain(vec_size);
        for (size_t k = 0; k < vec_size; k++)
        {
            Beta1_plain[k] = commRando.pop_bool();
        }
        if (myengine.party_id == 2)
        {
            myengine.privateCompare(Xbits_plain, Beta1_plain, R_plain, Beta2_plain);
        }
        else
        {
            myengine.privateCompare(X_bits[myengine.party_id], Beta1_plain, R_plain, Beta2_plain);
        }
        if (myengine.party_id == 2)
            for (size_t k = 0; k < vec_size; k++)
            {
                assert(((uint64_t)X_plain[k] > (uint64_t)R_plain[k]) == (Beta1_plain[k] ^ Beta2_plain[k]));
            }
        // printf("X_plain\n");
        // printVecI64(X_plain);
        // printf("R_plain\n");
        // printVecI64(R_plain);
        // printf("Beta_plain\n");
        // for (size_t k = 0; k < vec_size; k++)
        // {
        //     printf("%d ", Beta1_plain[k] ^ Beta2_plain[k]);
        // }
        // printf("\n");
    }
    {
        size_t vec_size = 10;
        VecI64 X_plain(vec_size), Xj(vec_size);
        vecbool A_plain(vec_size), Aj(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = commRando.pop_int64() / 2;
        }
        X_plain[0] = 0;
        X_plain[1] = -1l;
        X_plain[2] = -2l;
        X_plain[3] = (int64_t)-1 << 62;
        X_plain[4] = ((int64_t)1 << 62) - 1;

        myengine.share(X_plain, Xj);
        Xj[2] = -1l;
        myengine.msb(Xj, Aj);
        myengine.reveal_vec_bools(Aj, A_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(A_plain[i] == (X_plain[i] < 0));
        }

        VecI64 Zj(vec_size), Z_plain(vec_size);
        myengine.bit_inject(Xj, Aj, Zj);
        myengine.reveal(Zj, Z_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(Z_plain[i] == (X_plain[i] * A_plain[i]));
        }
        
        // printf("[LOG]:\tstart test4\n");
        // print_plus();
        // printf("Xj\n");
        // printVecI64(Xj);
        // printf("X_plain\n");
        // printVecI64(X_plain);
        // printf("A_plain\n");
        // printVecbool(A_plain);
        // printf("Z_plain\n");
        // printVecI64(Z_plain);
        // print_minus();
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
    MyEngine myengine(&socketMessenger, party_id);
    test1(myengine);
    test2(myengine);
    return 0;
}
