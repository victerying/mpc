#include "SecurennEngine.hpp"
#include "omp.h"
#include "bitset"
#include "../utils/util.hpp"

using std::bitset;

#define MYDEBUG

// #if defined MYDEBUG

// #endif

/// @brief 两个输入都在 [0, L-2]范围内
/// @param a
/// @param b
/// @return
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

SecurennEngine::SecurennEngine(VirtuleMessenger *_messenger, size_t _party_id) : messenger(_messenger), party_id(_party_id), prgs()
{
    this->message_send.resize(3);
    this->message_recv.resize(3);

    switch (this->party_id)
    {
    case 0:
        prgs.push_back(PRgenerator(KEY0, 0)); // PRG for self
        prgs.push_back(PRgenerator(KEY01, 0));
        prgs.push_back(PRgenerator(KEY20, 0));
        break;
    case 1:
        prgs.push_back(PRgenerator(KEY01, 0));
        prgs.push_back(PRgenerator(KEY1, 0)); // PRG for self
        prgs.push_back(PRgenerator(KEY12, 0));
        break;
    case 2:
        prgs.push_back(PRgenerator(KEY20, 0));
        prgs.push_back(PRgenerator(KEY12, 0));
        prgs.push_back(PRgenerator(KEY2, 0)); // PRG for self
        break;
    default:
        break;
    }
}

SecurennEngine::~SecurennEngine()
{
}

void SecurennEngine::matrix_mul(const MatrixXl &Xj, const MatrixXl &Yj, MatrixXl &Zj)
{
    size_t rows = Xj.rows();
    size_t lens = Xj.cols();
    assert(Yj.rows() == lens);
    size_t cols = Yj.cols();

    if (party_id == 2)
    {
        MatrixXl A[2];
        MatrixXl B[2];
        MatrixXl C[2];
        for (size_t i = 0; i < 2; i++)
        {
            A[i].resize(rows, lens);
            for (size_t j = 0; j < rows; j++)
            {
                for (size_t k = 0; k < lens; k++)
                {
                    A[i](j, k) = prgs[2].int_buff[i * rows * lens + j * lens + k];
                }
            }
        }
        int temp_index = 2 * rows * lens;
        for (size_t i = 0; i < 2; i++)
        {
            B[i].resize(lens, cols);
            for (size_t j = 0; j < lens; j++)
            {
                for (size_t k = 0; k < cols; k++)
                {
                    B[i](j, k) = prgs[2].int_buff[temp_index + i * lens * cols + j * cols + k];
                }
            }
        }
        temp_index += 2 * lens * cols;

        C[0].resize(rows, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t k = 0; k < cols; k++)
            {
                C[0](i, k) = prgs[2].int_buff[temp_index + i * cols + k];
            }
        }
        MatrixXl C_plain = (A[0] + A[1]) * (B[0] + B[1]);
        C[1] = C_plain - C[0];
        for (size_t i = 0; i < 2; i++)
        {
            matrix2string(A[i], message_send[i]);

            string temp_str;
            matrix2string(B[i], temp_str);
            message_send[i] = message_send[i] + temp_str;
            matrix2string(C[i], temp_str);
            message_send[i] = message_send[i] + temp_str;
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        message_send[0].clear();
        message_send[1].clear();
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信
        Zj.resize(rows, cols);
    }

    else // for j = 0,1  (party_id)
    {
        MatrixXl Aj(rows, lens);
        MatrixXl Bj(lens, cols);
        MatrixXl Cj(rows, cols);
        for (size_t i = 0; i < 3; i++)
        {
            message_send[i].clear();
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string str_a, str_b, str_c;
        size_t temp_index = 0;
        str_a = message_recv[2].substr(temp_index, rows * lens * sizeof(int64_t));

        temp_index += rows * lens * sizeof(int64_t);
        str_b = message_recv[2].substr(temp_index, lens * cols * sizeof(int64_t));

        temp_index += lens * cols * sizeof(int64_t);
        str_c = message_recv[2].substr(temp_index, rows * cols * sizeof(int64_t));
        assert(str_c.size() == rows * cols * sizeof(int64_t));

        string2matrix(str_a, Aj);
        string2matrix(str_b, Bj);
        string2matrix(str_c, Cj);

        MatrixXl Ej = Xj - Aj;
        MatrixXl Fj = Yj - Bj;
        string temp1, temp2;
        matrix2string(Ej, temp1);
        matrix2string(Fj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear();                              // 减少带宽
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信

        temp1 = message_recv[1 - party_id].substr(0, temp1.size());
        temp2 = message_recv[1 - party_id].substr(temp1.size(), temp2.size());
        MatrixXl E_plain(rows, lens);
        MatrixXl F_plain(lens, cols);
        string2matrix(temp1, E_plain);
        string2matrix(temp2, F_plain);
        E_plain = E_plain + Ej;
        F_plain = F_plain + Fj;
        Zj = Cj + (E_plain * Bj) + (Aj * F_plain) + (E_plain * F_plain * (int64_t)party_id);
    }
}

void SecurennEngine::vector_mul(const vector<int64_t> &Xj, const vector<int64_t> &Yj, vector<int64_t> &Zj)
{
    size_t vec_size = Xj.size();
    assert(Yj.size() == vec_size);
    Zj.resize(vec_size);
    if (party_id == 2)
    {
        vector<int64_t> A[2];
        vector<int64_t> B[2];
        vector<int64_t> C[2];
        for (size_t i = 0; i < 2; i++)
        {
            A[i].resize(vec_size);
            B[i].resize(vec_size);
            C[i].resize(vec_size);
        }
        for (size_t i = 0; i < 2; i++)
        {
            for (size_t j = 0; j < vec_size; j++)
            {
                A[i][j] = prgs[2].int_buff[i * vec_size * 3 + j * 3];
                B[i][j] = prgs[2].int_buff[i * vec_size * 3 + j * 3 + 1];
                C[i][j] = prgs[2].int_buff[i * vec_size * 3 + j * 3 + 2];
            }
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            C[1][i] = (A[0][i] + A[1][i]) * (B[0][i] + B[1][i]) - C[0][i];
        }

        string tempstr;
        for (size_t i = 0; i < 2; i++)
        {
            message_send[i].clear();
            vector2string(A[i], tempstr);
            message_send[i] += tempstr;
            vector2string(B[i], tempstr);
            message_send[i] += tempstr;
            vector2string(C[i], tempstr);
            message_send[i] += tempstr;
            assert(message_send[i].size() == 3 * vec_size * sizeof(int64_t));
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        message_send[0].clear();
        message_send[1].clear();
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信

#if defined MYDEBUG
        vector<int64_t> dummy(vec_size);
        vector<int64_t> A_plain, B_plain, C_plain, E_plain, F_plain, X_plain, Y_plain, Z_plain;
        reveal(dummy, A_plain);
        reveal(dummy, B_plain);
        reveal(dummy, C_plain);
        reveal(dummy, E_plain);
        reveal(dummy, F_plain);
        reveal(dummy, X_plain);
        reveal(dummy, Y_plain);
        reveal(dummy, Z_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(A_plain[i] == A[0][i] + A[1][i]);
            assert(B_plain[i] == B[0][i] + B[1][i]);

            assert(C_plain[i] == C[0][i] + C[1][i]);
            assert(A_plain[i] * B_plain[i] == C_plain[i]);
        }

#endif
    }

    else // for j = 0,1  (party_id)
    {
        vector<int64_t> Aj(vec_size);
        vector<int64_t> Bj(vec_size);
        vector<int64_t> Cj(vec_size);
        vector<int64_t> Ej(vec_size);
        vector<int64_t> Fj(vec_size);
        for (size_t i = 0; i < 3; i++)
        {
            message_send[i].clear();
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string str_a, str_b, str_c;
        size_t temp_index = 0;
        size_t str_size = vec_size * sizeof(int64_t);
        assert(message_recv[2].size() == 3 * vec_size * sizeof(int64_t));

        str_a = message_recv[2].substr(temp_index, str_size);
        temp_index += str_size;
        str_b = message_recv[2].substr(temp_index, str_size);
        temp_index += str_size;
        str_c = message_recv[2].substr(temp_index, str_size);
        assert(str_c.size() == vec_size * sizeof(int64_t));
        string2vector(str_a, Aj);
        string2vector(str_b, Bj);
        string2vector(str_c, Cj);

        for (size_t i = 0; i < vec_size; i++)
        {
            Ej[i] = Xj[i] - Aj[i];
            Fj[i] = Yj[i] - Bj[i];
        }
        string temp1, temp2;
        vector2string(Ej, temp1);
        vector2string(Fj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear();                              // 减少带宽
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信
        assert(message_recv[1 - party_id].size() == temp1.size() + temp2.size());

        temp1 = message_recv[1 - party_id].substr(0, temp1.size());
        temp2 = message_recv[1 - party_id].substr(temp1.size(), temp2.size());
        vector<int64_t> E_plain(vec_size);
        vector<int64_t> F_plain(vec_size);
        string2vector(temp1, E_plain);
        string2vector(temp2, F_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            E_plain[i] = Ej[i] + E_plain[i];
            F_plain[i] = Fj[i] + F_plain[i];
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            Zj[i] = Cj[i] + (E_plain[i] * Bj[i]) + (F_plain[i] * Aj[i]) + (F_plain[i] * E_plain[i] * (int64_t)party_id);
        }
#if defined MYDEBUG
        vector<int64_t> A_plain, B_plain, C_plain, debug_E_plain, debug_F_plain, X_plain, Y_plain, Z_plain;
        reveal(Aj, A_plain);
        reveal(Bj, B_plain);
        reveal(Cj, C_plain);
        reveal(Ej, debug_E_plain);
        reveal(Fj, debug_F_plain);
        reveal(Xj, X_plain);
        reveal(Yj, Y_plain);
        reveal(Zj, Z_plain);
        for (size_t i = 0; i < vec_size; i++)
        {

            assert(debug_E_plain[i] == X_plain[i] - A_plain[i]);
            assert(debug_F_plain[i] == Y_plain[i] - B_plain[i]);
            assert(E_plain[i] == X_plain[i] - A_plain[i]);
            assert(F_plain[i] == Y_plain[i] - B_plain[i]);
            assert(Z_plain[i] == X_plain[i] *Y_plain[i]);
        }
#endif
    }
}

void SecurennEngine::share(const vector<int64_t> &X_plain, vector<int64_t> &Xj)
{
    size_t vec_size = X_plain.size();
    Xj.resize(vec_size);
    if (party_id == 2)
    {
    }
    else
    {
        for (size_t i = 0; i < vec_size; i++)
        {
            int64_t temp = prgs[1 - party_id].pop_int64();
            Xj[i] = X_plain[i] * party_id;
            Xj[i] = party_id == 0 ? Xj[i] - temp : Xj[i] + temp;
        }
    }
}

void SecurennEngine::reveal(const vector<int64_t> &Xj, vector<int64_t> &X_plain)
{
    size_t vec_size = Xj.size();
    X_plain.resize(vec_size);
    if (party_id == 2)
    {
        for (size_t i = 0; i < 2; i++)
        {
            message_send[i].clear();
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        vector<int64_t> X0, X1;
        string2vector(message_recv[0], X0);
        string2vector(message_recv[1], X1);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = X0[i] + X1[i];
        }
    }
    else
    {

        vector2string(Xj, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2vector(message_recv[1 - party_id], X_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] += Xj[i];
        }
    }
}

void SecurennEngine::reveal_L_1(const vector<int64_t> &Xj, vector<int64_t> &X_plain)
{
    size_t vec_size = Xj.size();
    X_plain.resize(vec_size);
    if (party_id == 2)
    {
        for (size_t i = 0; i < 2; i++)
        {
            message_send[i].clear();
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        vector<int64_t> X0, X1;
        string2vector(message_recv[0], X0);
        string2vector(message_recv[1], X1);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = add_L_1(X0[i], X1[i]);
        }
    }
    else
    {
        vector2string(Xj, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2vector(message_recv[1 - party_id], X_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = add_L_1(Xj[i], X_plain[i]);
        }
    }
}

void SecurennEngine::reveal_zp(const vector<zpshare> &Xbits_j, vector<zpshare> &Xbits_plain)
{
    size_t vec_size = Xbits_j.size();
    Xbits_plain.resize(vec_size);
    if (party_id == 2)
    {
        for (size_t i = 0; i < 2; i++)
        {
            message_send[i].clear();
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        vector<zpshare> Xbits_0, Xbits_1;
        string2vecZpshare(message_recv[0], Xbits_0);
        string2vecZpshare(message_recv[1], Xbits_1);
        for (size_t i = 0; i < vec_size; i++)
        {
            for (size_t k = 0; k < 64; k++)
            {
                Xbits_plain[i][k] = (Xbits_0[i][k] + Xbits_1[i][k]) % P;
            }
        }
    }
    else
    {
        vecZpshare2string(Xbits_j, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2vecZpshare(message_recv[1 - party_id], Xbits_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            for (size_t k = 0; k < 64; k++)
            {
                Xbits_plain[i][k] = (Xbits_j[i][k] + Xbits_plain[i][k]) % P;
            }
        }
    }
}

/// @brief Securenn 的PrivateCompare协议
/// @param Xbits_j zpshare[i] 表示x从小到达第i个比特在Zp上的秘密分享，party2秩序给出一段填充的数据
/// @param Beta1_plain 明文，party0,1输入必须相同，party2秩序给出一段填充的数据
/// @param R_palin 明文，party0,1输入必须相同，party2秩序给出一段填充的数据
/// @param Beta2_plain 输出，明文，被party2知道
void SecurennEngine::privateCompare(const vector<zpshare> &Xbits_j, const vector<int8_t> &Beta1_plain, const vector<int64_t> &R_palin, vector<int8_t> &Beta2_plain)
{
    size_t vec_size = Xbits_j.size();
    assert(Beta1_plain.size() == vec_size);
    assert(R_palin.size() == vec_size);
    Beta2_plain.resize(vec_size);
    if (party_id == 0 || party_id == 1)
    {
        vector<zpshare> S_plain, U_plain;
        S_plain.resize(vec_size);
        U_plain.resize(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            for (size_t j = 0; j < 64; j++)
            {
                S_plain[i][j] = prgs[1 - party_id].int_buff[i * 64 + j];
                U_plain[i][j] = prgs[1 - party_id].int_buff[vec_size * 64 + i * 64 + j];
                S_plain[i][j] = S_plain[i][j] % P;
                if (0 == S_plain[i][j])
                {
                    S_plain[i][j] = 1;
                }
                U_plain[i][j] = U_plain[i][j] % P;
                if (0 == U_plain[i][j])
                {
                    U_plain[i][j] = 1;
                }
            }
        }
        vector<zpshare> Dbits_j(vec_size);
        vector<zpshare> Rbits_plain(vec_size);
        vector<zpshare> Tbits_plain(vec_size);
        vector<zpshare> Wbits_j(vec_size);
        vector<zpshare> Cbits_j(vec_size);

        for (size_t i = 0; i < vec_size; i++)
        {

            bitset<64> bstr(R_palin[i]);
            bitset<64> bstt(R_palin[i] + 1);
            for (size_t k = 0; k < 64; k++)
            {
                Rbits_plain[i][k] = (bool)bstr[k];
                Tbits_plain[i][k] = (bool)bstt[k];
            }
            uint16_t sum_w = 0;
            for (size_t k = 63; k != -1; k--)
            {
                if (Beta1_plain[i] == 0)
                {
                    // r[k] + 1 - x
                    Cbits_j[i][k] = (Xbits_j[i][k] * (P - 1)) % P + (uint16_t)party_id * Rbits_plain[i][k] + (uint16_t)party_id + sum_w;
                    Cbits_j[i][k] = Cbits_j[i][k] % P;
                    Wbits_j[i][k] = Xbits_j[i][k] + (uint16_t)party_id * Rbits_plain[i][k] + (P - 2) * Rbits_plain[i][k] * Xbits_j[i][k];
                    Wbits_j[i][k] = Wbits_j[i][k] % P;
                    sum_w = sum_w + Wbits_j[i][k];
                }
                else if (Beta1_plain[i] == 1 && R_palin[i] != -1)
                {
                    Cbits_j[i][k] = (P - 1) * party_id * Tbits_plain[i][k] + Xbits_j[i][k] + (uint16_t)party_id + sum_w;
                    Cbits_j[i][k] = Cbits_j[i][k] % P;
                    Wbits_j[i][k] = Xbits_j[i][k] + (uint16_t)party_id * Tbits_plain[i][k] + (P - 2) * Tbits_plain[i][k] * Xbits_j[i][k];
                    Wbits_j[i][k] = Wbits_j[i][k] % P;
                    sum_w = sum_w + Wbits_j[i][k];
                }
                else
                {
                    printf("track 3. ");
                    if (k != 0)
                    {
                        Cbits_j[i][k] = (1 - party_id) * (1 + U_plain[i][k]) + (P - party_id) * U_plain[i][k];
                        Cbits_j[i][k] = Cbits_j[i][k] % P;
                    }
                    else
                    {
                        Cbits_j[i][k] = party_id == 0 ? U_plain[i][k] : P - U_plain[i][k];
                    }
                }
            }

            for (size_t k = 0; k < 64; k++)
            {
                Dbits_j[i][k] = S_plain[i][k] * Cbits_j[i][k];
                Dbits_j[i][k] = Dbits_j[i][k] % P;
            }
        }
        // 这里应该置换Dbits_j
        message_send[1 - party_id].clear();
        vecZpshare2string(Dbits_j, message_send[2]);
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
#if defined MYDEBUG
        // printf("privateCompare debuging\n");
        // printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        vector<zpshare> Xbits_debug(vec_size);
        vector<zpshare> Wbits_debug(vec_size);
        vector<zpshare> Cbits_debug(vec_size);
        reveal_zp(Xbits_j, Xbits_debug);
        reveal_zp(Wbits_j, Wbits_debug);
        reveal_zp(Cbits_j, Cbits_debug);
        

        // printf("Xbits_debug:\n");
        // print_vec_zpshare(Xbits_debug);
        // printf("Rbits_plain:\n");
        // print_vec_zpshare(Rbits_plain);
        // printf("Tbits_plain:\n");
        // print_vec_zpshare(Tbits_plain);
        // printf("Wbits_debug:\n");
        // print_vec_zpshare(Wbits_debug);
        // printf("Cbits_debug:\n");
        // print_vec_zpshare(Cbits_debug);

        for (size_t i = 0; i < vec_size; i++)
        {
            uint16_t sum_w = 0;
            for (size_t k = 63; k != -1; k--)
            {
                if (Beta1_plain[i] == 0)
                {
                    assert(Wbits_debug[i][k] == Xbits_debug[i][k] ^ Rbits_plain[i][k]);
                    assert(Cbits_debug[i][k] == Rbits_plain[i][k] + 1 - Xbits_debug[i][k] + sum_w);
                    sum_w = sum_w + Wbits_debug[i][k];
                }
                else if (Beta1_plain[i] == 1 && R_palin[i] != -1)
                {
                    assert(Wbits_debug[i][k] == Xbits_debug[i][k] ^ Tbits_plain[i][k]);
                    assert(Cbits_debug[i][k] == Xbits_debug[i][k] + 1 - Tbits_plain[i][k] + sum_w);
                    sum_w = sum_w + Wbits_debug[i][k];
                }
                else
                {
                    assert(Cbits_debug[i][k] == (k == 0 ? 0 : 1));
                }
            }
        }

#endif
    }
    else // party_id == 2
    {
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        vector<zpshare> Dbits_0, Dbits_1;
        string2vecZpshare(message_recv[0], Dbits_0);
        string2vecZpshare(message_recv[1], Dbits_1);
        assert(Dbits_0.size() == vec_size);
        assert(Dbits_1.size() == vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            zpshare dbits_plain;
            Beta2_plain[i] = 0;
            for (size_t k = 0; k < 64; k++)
            {
                dbits_plain[k] = Dbits_0[i][k] + Dbits_1[i][k];
                dbits_plain[k] = dbits_plain[k] % P;
                if (0 == dbits_plain[k])
                {
                    Beta2_plain[i] = 1;
                    break;
                }
            }
        }
#if defined MYDEBUG
        vector<zpshare> Xbits_debug(vec_size);
        vector<zpshare> Wbits_debug(vec_size);
        vector<zpshare> Cbits_debug(vec_size);
        vector<zpshare> dummy_pad(vec_size);
        reveal_zp(dummy_pad, Xbits_debug);
        reveal_zp(dummy_pad, Wbits_debug);
        reveal_zp(dummy_pad, Cbits_debug);
#endif
    }
}

/// @brief share convert
/// @param Aj Z_L上的秘密分享
/// @param Yj Z_L_1上的秘密分享
void SecurennEngine::shareConvert(const vector<int64_t> &Aj, vector<int64_t> &Yj)
{
    size_t vec_size = Aj.size();
    Yj.resize(vec_size);
    if (2 != party_id)
    {
        // 公共随机数
        vector<int8_t> Eta2_plain(vec_size);
        vector<int64_t> R_plain(vec_size);
        vector<int64_t> R_share[2];
        R_share[0].resize(vec_size);
        R_share[1].resize(vec_size);
        vector<int8_t> Alpha_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Eta2_plain[i] = prgs[1 - party_id].bool_buff[i];
            R_share[0][i] = prgs[1 - party_id].int_buff[i];
            R_share[1][i] = prgs[1 - party_id].int_buff[vec_size + i];
            R_plain[i] = R_share[0][i] + R_share[1][i];
            Alpha_plain[i] = (uint64_t)R_plain[i] < (uint64_t)R_share[0][i];
        }
        vector<int64_t> Xj(vec_size);

        vector<int8_t> Beta_j(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            // x = a + r
            Xj[i] = Aj[i] + R_share[party_id][i];
            Beta_j[i] = (uint64_t)Xj[i] < (uint64_t)R_share[party_id][i];
        }
        message_send[1 - party_id].clear();
        vector2string(Xj, message_send[2]);
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        message_send[2].clear();
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信

        vector<zpshare> Xbits_j;
        // 在 Z_L_1上的秘密分享
        vector<int64_t> Delta_j;
        size_t temp_index = vec_size * 64 * sizeof(uint16_t);
        string2vecZpshare(message_recv[2].substr(0, temp_index), Xbits_j);
        string2vector(message_recv[2].substr(temp_index, message_recv[2].size() - temp_index), Delta_j);
        vector<int8_t> dummy_beta2plain;
        privateCompare(Xbits_j, Eta2_plain, R_plain, dummy_beta2plain); // 第三轮通信
        message_send[1 - party_id].clear();
        message_send[2].clear();
        messenger->send_and_recv(message_send, message_recv); // 第四轮通信
        vector<int64_t> Eta1_j;
        // 在 Z_L_1上的秘密分享
        string2vector(message_recv[2], Eta1_j);
        vector<int64_t> Eta_j(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Eta_j[i] = add_L_1(Eta1_j[i], party_id * Eta2_plain[i]);
            uint64_t temp = add_L_1(Eta1_j[i], Eta1_j[i]);
            temp = temp * Eta2_plain[i];
            Eta_j[i] = sub_L_1(Eta_j[i], temp);
        }
        vector<int64_t> Theta_j(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            uint64_t tempj = (1 - party_id) * (L_1 - 1 - Alpha_plain[i]);
            Theta_j[i] = add_L_1(Beta_j[i], tempj);
            Theta_j[i] = add_L_1(Theta_j[i], Eta_j[i]);
            Theta_j[i] = add_L_1(Theta_j[i], Delta_j[i]);
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            Yj[i] = sub_L_1(Aj[i], Theta_j[i]);
        }
#if defined MYDEBUG
        vector<int64_t> A_debug, X_debug;
        vector<zpshare> Xbits_debug;
        vector<int64_t> Delta_debug;
        vector<int64_t> Eta__debug;
        vector<int64_t> Theta_debug;
        reveal(Aj, A_debug);
        reveal(Xj, X_debug);
        reveal_zp(Xbits_j, Xbits_debug);
        reveal_L_1(Delta_j, Delta_debug);
        reveal_L_1(Eta_j, Eta__debug);
        reveal_L_1(Theta_j, Theta_debug);
        // printf("Aj:\n");
        // print_vec_int64(Aj);
        // printf("R_share[%lu]\n", party_id);
        // print_vec_int64(R_share[party_id]);
        // printf("Xj:\n");
        // print_vec_int64(Xj);

        // printf("A_debug:\n");
        // print_vec_int64(A_debug);
        // printf("R_plain:\n");
        // print_vec_int64(R_plain);
        // printf("X_debug:\n");
        // print_vec_int64(X_debug);

        // printf("Alpha_plain:\n");
        // for (size_t i = 0; i < vec_size; i++)
        // {
        //     printf("%d ", Alpha_plain[i]);
        // }
        // printf("\n");

        // printf("Beta_j:\n");
        // for (size_t i = 0; i < vec_size; i++)
        // {
        //     printf("%d ", Beta_j[i]);
        // }
        // printf("\n");

        // printf("Delta_debug:\n");
        // print_vec_int64(Delta_debug);
        // printf("Eta__debug:\n");
        // print_vec_int64(Eta__debug);
        // printf("Theta_debug:\n");
        // print_vec_int64(Theta_debug);

#endif
    }
    else // 2 == party_id
    {
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        vector<int64_t> X0, X1, X_plain;
        X_plain.resize(vec_size);
        vector<int8_t> Delta_plain(vec_size);
        string2vector(message_recv[0], X0);
        string2vector(message_recv[1], X1);
        assert(X0.size() == vec_size && X1.size() == vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = X0[i] + X1[i];
            Delta_plain[i] = (uint64_t)X_plain[i] < (uint64_t)X0[i];
        }
        vector<zpshare> Xbits_plain(vec_size), Xbits_0(vec_size), Xbits_1(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            int64_t &x_plain = X_plain[i];
            zpshare &xbits_plain = Xbits_plain[i];
            zpshare &xbits_0 = Xbits_0[i];
            zpshare &xbits_1 = Xbits_1[i];
            bitset<64> btsx(x_plain);
            for (size_t k = 0; k < 64; k++)
            {
                xbits_plain[k] = (bool)btsx[k];
                xbits_0[k] = this->prgs[2].int_buff[i * 64 + k];
                xbits_0[k] = xbits_0[k] % P;
                xbits_1[k] = xbits_plain[k] + P - xbits_0[k];
                xbits_1[k] = xbits_1[k] % P;
            }
        }
        // 在 Z_L_1上的秘密分享
        vector<int64_t> Delta_0(vec_size), Delta_1(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Delta_0[i] = this->prgs[2].int_buff[i];
            Delta_0[i] = (uint64_t)Delta_0[i] % L_1;
            uint64_t temp = Delta_plain[i];
            Delta_1[i] = sub_L_1(temp, Delta_0[i]);
        }
        string temp_str;
        vecZpshare2string(Xbits_0, message_send[0]);
        vector2string(Delta_0, temp_str);
        message_send[0] = message_send[0] + temp_str;
        vecZpshare2string(Xbits_1, message_send[1]);
        vector2string(Delta_1, temp_str);
        message_send[1] = message_send[1] + temp_str;
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信
        vector<int8_t> Eta1_plain;
        // 前三个输入只是为了占位
        privateCompare(Xbits_plain, Delta_plain, X_plain, Eta1_plain); // 第三轮通信
        // 在 Z_L_1上的秘密分享
        vector<int64_t> Eta1_0(vec_size), Eta1_1(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Eta1_0[i] = this->prgs[2].int_buff[vec_size + i];
            Eta1_0[i] = Eta1_0[i] % L_1;
            uint64_t temp = Eta1_plain[i];
            Eta1_1[i] = sub_L_1(temp, Eta1_0[i]);
        }
        vector2string(Eta1_0, message_send[0]);
        vector2string(Eta1_1, message_send[1]);
        messenger->send_and_recv(message_send, message_recv); // 第四轮通信
#if defined MYDEBUG
        vector<int64_t> dummy_pad(vec_size), dummy_pad1;
        vector<int64_t> A_debug, X_debug;
        vector<zpshare> Xbits_debug;
        vector<int64_t> Theta_debug;
        reveal(Aj, A_debug);
        reveal(dummy_pad, X_debug);
        reveal_zp(Xbits_plain, Xbits_debug);
        reveal_L_1(dummy_pad, dummy_pad1);
        reveal_L_1(dummy_pad, dummy_pad1);
        reveal_L_1(dummy_pad, Theta_debug);

        // printf("Delta_plain:\n");
        // for (size_t i = 0; i < vec_size; i++)
        // {
        //     printf("%d ", Delta_plain[i]);
        // }
        // printf("\n");  

#endif
    }
}

/// @brief 计算msb
/// @param Aj Z_L_1上的秘密分享
/// @param Bj Z_L上的秘密分享
void SecurennEngine::msb(const vector<int64_t> &Aj, vector<int64_t> &Bj)
{
    size_t vec_size = Aj.size();
    Bj.resize(vec_size);
    if (party_id == 2)
    {
        // Z_L_1上的秘密分享
        vector<int64_t> X_plain(vec_size), X_0(vec_size), X_1(vec_size);
        // Z_L上的秘密分享
        vector<int64_t> Xlsb_0(vec_size), Xlsb_1(vec_size);
        vector<zpshare> Xbits_0(vec_size), Xbits_1(vec_size), Xbits_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_0[i] = this->prgs[2].int_buff[i];
            X_0[i] = X_0[i] % L_1;
            X_1[i] = this->prgs[2].int_buff[vec_size + i];
            X_1[i] = X_1[i] % L_1;
            X_plain[i] = add_L_1(X_0[i], X_1[i]);
            bitset<64> btst_x(X_plain[i]);
            Xlsb_0[i] = this->prgs[2].int_buff[2 * vec_size + i];
            Xlsb_1[i] = (bool)btst_x[0] - Xlsb_0[i];
            for (size_t k = 0; k < 64; k++)
            {
                Xbits_plain[i][k] = (bool)btst_x[k];
                Xbits_1[i][k] = this->prgs[2].int_buff[3 * vec_size + i * 64 + k];
                Xbits_1[i][k] = Xbits_1[i][k] % P;

                Xbits_0[i][k] = Xbits_plain[i][k] + P - Xbits_1[i][k];
                Xbits_0[i][k] = Xbits_0[i][k] % P;
            }
        }
        string temp_str;
        vector2string(X_0, message_send[0]);
        vector2string(Xlsb_0, temp_str);
        message_send[0] = message_send[0] + temp_str;
        vecZpshare2string(Xbits_0, temp_str);
        message_send[0] = message_send[0] + temp_str;

        vector2string(X_1, message_send[1]);
        vector2string(Xlsb_1, temp_str);
        message_send[1] = message_send[1] + temp_str;
        vecZpshare2string(Xbits_1, temp_str);
        message_send[1] = message_send[1] + temp_str;
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信

        message_send[0].clear();
        message_send[1].clear();
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信
        vector<int8_t> Beta1_plain, dummy_Beta0(vec_size);
        privateCompare(Xbits_plain, dummy_Beta0, X_plain, Beta1_plain); // 第三轮通信

        // Z_L上的秘密分享
        vector<int64_t> Beta1_0(vec_size), Beta1_1(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Beta1_0[i] = prgs[2].int_buff[4 * vec_size + i];
            Beta1_1[i] = Beta1_plain[i] - Beta1_0[i];
        }
        vector2string(Beta1_0, message_send[0]);
        vector2string(Beta1_1, message_send[1]);
        messenger->send_and_recv(message_send, message_recv); // 第四轮通信

        vector_mul(X_plain, X_plain, X_plain); // 第五轮通信
#if defined MYDEBUG
        vector<int64_t> dummy_pad(vec_size);
        //  加上R，这些都是在L_1上的秘密分享 Y = A + A R = Y + X   msb(A) == lsb(Y)
        vector<int64_t> A_debug, X_debug, Y_debug;
        // Beta_debug = X > R;   Delta_debug = lsb(X) xor lsb(R)  在L上的秘密分享
        vector<int64_t> Beta_debug, Delta_debug, Xlsb_debug, Theta_debug;

        reveal_L_1(Aj, A_debug);
        reveal_L_1(dummy_pad, X_debug);
        reveal_L_1(dummy_pad, Y_debug);

        reveal_L_1(dummy_pad, Beta_debug);
        reveal_L_1(dummy_pad, Delta_debug);
        reveal_L_1(dummy_pad, Xlsb_debug);
        reveal(dummy_pad, Theta_debug);

#endif
    }
    else // party_id == 0,1
    {
        // 公共随机数
        vector<int8_t> Beta0_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Beta0_plain[i] = prgs[1 - party_id].bool_buff[i];
        }

        message_send[2].clear();
        message_send[1 - party_id].clear();
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        // Xj: Z_L_1上的秘密分享, Xlsb_j:Z_L上的秘密分享
        vector<int64_t> Xj, Xlsb_j;
        vector<zpshare> Xbits_j;
        size_t temp_index = vec_size * sizeof(int64_t);
        string2vector(message_recv[2].substr(0, temp_index), Xj);
        string2vector(message_recv[2].substr(temp_index, temp_index), Xlsb_j);
        string2vecZpshare(message_recv[2].substr(2 * temp_index, message_recv[2].size() - temp_index), Xbits_j);

        // Z_L_1上的秘密分享
        vector<int64_t> Yj(vec_size);
        vector<int64_t> Rj(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Yj[i] = add_L_1(Aj[i], Aj[i]);
            Rj[i] = add_L_1(Yj[i], Xj[i]);
        }
        vector2string(Rj, message_send[1 - party_id]);
        messenger->send_and_recv(message_send, message_recv); // 第二轮通信
        vector<int64_t> R_plain;
        string2vector(message_recv[1 - party_id], R_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            R_plain[i] = add_L_1(R_plain[i], Rj[i]);
        }
        vector<int8_t> dummy_Beta1;
        privateCompare(Xbits_j, Beta0_plain, R_plain, dummy_Beta1); // 第三轮通信
        message_send[2].clear();
        message_send[1 - party_id].clear();
        messenger->send_and_recv(message_send, message_recv); // 第四轮通信
        vector<int64_t> Beta1_j;
        string2vector(message_recv[2], Beta1_j);
        vector<int64_t> R0_plain(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            R0_plain[i] = (uint64_t)R_plain[i] % 2;
        }
        vector<int64_t> Beta_j(vec_size), Delta_j(vec_size), Theta_j(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Beta_j[i] = Beta1_j[i] + party_id * Beta0_plain[i] - 2 * Beta0_plain[i] * Beta1_j[i];
            Delta_j[i] = Xlsb_j[i] + party_id * R0_plain[i] - 2 * R0_plain[i] * Xlsb_j[i];
        }
        vector_mul(Beta_j, Delta_j, Theta_j); // 第五轮通信
        for (size_t i = 0; i < vec_size; i++)
        {
            Bj[i] = Beta_j[i] + Delta_j[i] - 2 * Theta_j[i];
        }
#if defined MYDEBUG
        //  加上R，这些都是在L_1上的秘密分享 Y = A + A R = Y + X   msb(A) == lsb(Y)
        vector<int64_t> A_debug, X_debug, Y_debug;
        // Beta_debug = X > R;   Delta_debug = lsb(X) xor lsb(R)  Theta_debug = Beta_debug * Delta_debug在L上的秘密分享
        vector<int64_t> Beta_debug, Delta_debug, Xlsb_debug, Theta_debug;

        reveal_L_1(Aj, A_debug);
        reveal_L_1(Xj, X_debug);
        reveal_L_1(Yj, Y_debug);

        reveal(Beta_j, Beta_debug);
        reveal(Delta_j, Delta_debug);
        reveal(Xlsb_j, Xlsb_debug);
        reveal(Theta_j, Theta_debug);
        // printf("in msb++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===\n");
        // printf("A_debug:\n");
        // print_vec_int64(A_debug);
        // printf("X_debug:\n");
        // print_vec_int64(X_debug);
        // printf("Y_debug:\n");
        // print_vec_int64(Y_debug);
        // printf("R_plain:\n");
        // print_vec_int64(R_plain);

        
        // printf("Xlsb_debug:\n");
        // print_vec_int64(Xlsb_debug);
        // printf("R0_plain:\n");
        // print_vec_int64(R0_plain);
        // printf("Delta_debug:\n");
        // print_vec_int64(Delta_debug);
        // printf("Beta_debug:\n");
        // print_vec_int64(Beta_debug);
        // printf("Theta_debug:\n");
        // print_vec_int64(Theta_debug);
#endif
    }
}
