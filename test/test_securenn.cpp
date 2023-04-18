#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <bitset>
#include "../mpc_engine/SecurennEngine.hpp"
#include "../utils/PRgenerator.hpp"
#include "../utils/SocketMessenger.hpp"
#include "../utils/util.hpp"
using std::bitset;
using std::cout;

/// @brief 测试公共的为随机数生成器
/// @param securenn_eigine
void test1(SecurennEngine &securenn_eigine)
{
    size_t party_id = securenn_eigine.party_id;
    printf("party_id: %lu\n", party_id);
    if (party_id == 2)
    {
        for (size_t j = 0; j < 2; j++)
        {
            printf("comman rand(%lu,%lu):\n", j, party_id);
            for (size_t k = 0; k < 5; k++)
            {
                printf("%016lX ", securenn_eigine.prgs[j].pop_int64());
            }
            printf("\n");
        }
    }
    else
    {
        printf("comman rand(%lu,%lu):\n", party_id, 1 - party_id);
        for (size_t k = 0; k < 5; k++)
        {
            printf("%016lX ", securenn_eigine.prgs[1 - party_id].pop_int64());
        }
        printf("\n");

        printf("comman rand(%lu,%lu):\n", party_id, (uint64_t)2);
        for (size_t k = 0; k < 5; k++)
        {
            printf("%016lX ", securenn_eigine.prgs[2].pop_int64());
        }
        printf("\n");
    }
}

/// @brief 测试分享和揭秘
/// @param securenn_eigine
void test2(SecurennEngine &securenn_eigine)
{
    printf("start test2\n");
    size_t party_id = securenn_eigine.party_id;
    size_t vec_size = 10;
    vector<int64_t> X_plain(vec_size);
    vector<int64_t> Xj, _X_plain;
    for (size_t l = 0; l < 1; l++)
    {
        printf("starting loop %lu\n", l);
        if (party_id == 2)
        {
            securenn_eigine.share(X_plain, Xj);
            securenn_eigine.reveal(Xj, _X_plain);
            printf("finish reveal, reveal :\n");
            for (size_t i = 0; i < vec_size; i++)
            {
                printf("%016lX ", _X_plain[i]);
            }
            printf("\n");
        }
        else
        {
            printf("plain_text:\n");
            for (size_t i = 0; i < vec_size; i++)
            {
                X_plain[i] = securenn_eigine.prgs[1 - party_id].int_buff[i];
                printf("%016lX ", X_plain[i]);
            }
            printf("\n");

            securenn_eigine.share(X_plain, Xj);

            printf("finish share, share :\n");
            for (size_t i = 0; i < vec_size; i++)
            {
                printf("%016lX ", Xj[i]);
            }
            printf("\n");

            securenn_eigine.reveal(Xj, _X_plain);
            printf("finish reveal, reveal :\n");
            for (size_t i = 0; i < vec_size; i++)
            {
                printf("%016lX ", _X_plain[i]);
            }
            printf("\n");
        }
    }
    printf("end test2\n");
}

/// @brief 测试向量乘法和矩阵乘法
/// @param securenn_eigine
void test3(SecurennEngine &securenn_eigine)
{
    size_t party_id = securenn_eigine.party_id;
    // 测试vector_mul
    {
        size_t vec_size = 5;
        vector<int64_t> Xj(vec_size), Yj(vec_size), Zj;
        for (size_t i = 0; i < vec_size; i++)
        {
            Xj[i] = securenn_eigine.prgs[party_id].pop_int64();
            Yj[i] = securenn_eigine.prgs[party_id].pop_int64();
        }
        securenn_eigine.vector_mul(Xj, Yj, Zj);
        vector<int64_t> X_plain, Y_plain, Z_plain;
        securenn_eigine.reveal(Xj, X_plain);
        securenn_eigine.reveal(Yj, Y_plain);
        securenn_eigine.reveal(Zj, Z_plain);

        printf("x :\n");
        print_vec_int64(X_plain);
        printf("y :\n");
        print_vec_int64(Y_plain);

        printf("x*y :\n");
        for (size_t i = 0; i < vec_size; i++)
        {
            printf("%016lX ", X_plain[i] * Y_plain[i]);
        }
        printf("\n");

        printf("z :\n");
        print_vec_int64(Z_plain);
    }

    // 测试matrix_mul
    {
        size_t rows = 3;
        size_t lens = 4;
        size_t cols = 5;
        MatrixXl Xj(rows, lens), Yj(lens, cols), Zj;
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < lens; j++)
            {
                Xj(i, j) = securenn_eigine.prgs[party_id].pop_int64();
            }
        }
        for (size_t i = 0; i < lens; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                Yj(i, j) = securenn_eigine.prgs[party_id].pop_int64();
            }
        }
        securenn_eigine.matrix_mul(Xj, Yj, Zj);

        MatrixXl X_plain(rows, lens), Y_plain(lens, cols), Z_plain(rows, cols);
        vector<int64_t> temp_vecj, temp_vec_plain;

        matrix2vector(Xj, temp_vecj);
        securenn_eigine.reveal(temp_vecj, temp_vec_plain);
        vector2matrix(temp_vec_plain, X_plain);

        matrix2vector(Yj, temp_vecj);
        securenn_eigine.reveal(temp_vecj, temp_vec_plain);
        vector2matrix(temp_vec_plain, Y_plain);

        matrix2vector(Xj, temp_vecj);
        securenn_eigine.reveal(temp_vecj, temp_vec_plain);
        vector2matrix(temp_vec_plain, X_plain);

        matrix2vector(Zj, temp_vecj);
        securenn_eigine.reveal(temp_vecj, temp_vec_plain);
        vector2matrix(temp_vec_plain, Z_plain);

        printf("X_plain:\n");
        cout << X_plain << "\n";
        printf("Y_plain:\n");
        cout << Y_plain << "\n";
        printf("X_plain* Y_plain:\n");
        cout << X_plain * Y_plain << "\n";
        printf("Z_plain:\n");
        cout << Z_plain << "\n";
    }
}

/// @brief 测试privateCompare
/// @param securenn_eigine
void test4(SecurennEngine &securenn_eigine)
{
    printf("start test4+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    size_t party_id = securenn_eigine.party_id;
    size_t vec_size = 5;

    vector<int64_t> X_plain(vec_size);
    vector<int64_t> R_plain(vec_size);
    vector<int8_t> Beta1_plain(vec_size);
    vector<int8_t> Beta2_plain(vec_size);
    vector<zpshare> Xbits_plain(vec_size);
    vector<zpshare> Xbits_share[2];
    Xbits_share[0].resize(vec_size);
    Xbits_share[1].resize(vec_size);
    if (2 != party_id)
    {
        size_t rerandom = 14;
        for (size_t i = 0; i < rerandom; i++)
        {
            securenn_eigine.prgs[1 - party_id].pop_int64();
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = securenn_eigine.prgs[1 - party_id].pop_int64();
            R_plain[i] = securenn_eigine.prgs[1 - party_id].pop_int64();
            if (i % 5 == 2)
            {
                R_plain[i] = -1;
            }
            Beta1_plain[i] = securenn_eigine.prgs[1 - party_id].pop_bool();
            bitset<64> btsx(X_plain[i]);
            for (size_t k = 0; k < 64; k++)
            {
                Xbits_plain[i][k] = (bool)btsx[k];
                Xbits_share[0][i][k] = securenn_eigine.prgs[1 - party_id].pop_int64();
                Xbits_share[0][i][k] = Xbits_share[0][i][k] % P;
                Xbits_share[1][i][k] = Xbits_plain[i][k] + P - Xbits_share[0][i][k];
                Xbits_share[1][i][k] = Xbits_share[1][i][k] % P;
            }
        }
        vector<zpshare> _Xbits_plain;
        securenn_eigine.reveal_zp(Xbits_share[party_id], _Xbits_plain);
        assert(Xbits_plain == _Xbits_plain);

        securenn_eigine.privateCompare(Xbits_share[party_id], Beta1_plain, R_plain, Beta2_plain);
        printf("X_plain:\n");
        print_vec_int64(X_plain);
        printf("R_plain:\n");
        print_vec_int64(R_plain);
        printf("Xbits_plain:\n");
        print_vec_zpshare(Xbits_plain);
        printf("Xbits_share[%lu]:\n", party_id);
        print_vec_zpshare(Xbits_share[party_id]);
        printf("Beta1_plain:\n");
        for (size_t i = 0; i < Beta1_plain.size(); i++)
        {
            printf("%d ", Beta1_plain[i]);
        }
        printf("\n");
        printf("X > R:\n");
        for (size_t i = 0; i < vec_size; i++)
        {
            printf("%d ", (uint64_t)X_plain[i] > (uint64_t)R_plain[i]);
        }
        printf("\n");
    }
    else
    {
        vector<zpshare> _Xbits_plain;
        securenn_eigine.reveal_zp(Xbits_plain, _Xbits_plain);

        securenn_eigine.privateCompare(Xbits_plain, Beta1_plain, R_plain, Beta2_plain);
        printf("Beta2_plain:\n");
        for (size_t i = 0; i < Beta2_plain.size(); i++)
        {
            printf("%d ", Beta2_plain[i]);
        }
        printf("\n");
    }
}

/// @brief 测试shareConvert和msb
/// @param securenn_eigine
void test5(SecurennEngine &securenn_eigine)
{
    printf("start test5+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    size_t party_id = securenn_eigine.party_id;
    size_t vec_size = 20;
    vector<int64_t> A_plain(vec_size);
    vector<int64_t> Aj, Yj, Y_plain;
    
    

    if (party_id == 2)
    {
        securenn_eigine.share(A_plain, Aj);
        securenn_eigine.shareConvert(Aj, Yj);
        securenn_eigine.reveal_L_1(Yj, Y_plain);
    }
    else
    {
        size_t rerandom = 14;
        for (size_t i = 0; i < rerandom; i++)
        {
            securenn_eigine.prgs[1 - party_id].pop_int64();
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            A_plain[i] = securenn_eigine.prgs[1 - party_id].pop_int64();
            if (A_plain[i] == -1)
            {
                A_plain[i] = 1;
            }
            
        }
        securenn_eigine.share(A_plain, Aj);
        securenn_eigine.shareConvert(Aj, Yj);
        securenn_eigine.reveal_L_1(Yj, Y_plain);
    }
    printf("A_plain\n");
    print_vec_int64(A_plain);
    printf("Y_plain\n");
    print_vec_int64(Y_plain);
    if (2 != party_id)
    {
        assert(A_plain == Y_plain);
    }
    
    vector<int64_t> Bj, B_plain;
    securenn_eigine.msb(Yj, Bj);
    securenn_eigine.reveal(Bj, B_plain);
    printf("B_plain\n");
    print_vec_int64(B_plain);

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
    SecurennEngine securenn_eigine(&socketMessenger, party_id);
    printf("start test securenn\n");
    test1(securenn_eigine);
    test2(securenn_eigine);
    test3(securenn_eigine);
    test4(securenn_eigine);
    test5(securenn_eigine);
    return 0;
}
