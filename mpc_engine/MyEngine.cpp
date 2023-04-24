#include "MyEngine.hpp"

// #define MYDEBUG

MyEngine::MyEngine(VirtuleMessenger *_messenger, size_t _party_id) : messenger(_messenger), party_id(_party_id)
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

MyEngine::~MyEngine()
{
}

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

void MyEngine::share(const VecI64 &X_plain, VecI64 &Xj)
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

void MyEngine::share_matrix(const MatrixXl &X_plain, MatrixXl &Xj)
{
    size_t rows = X_plain.rows(), cols = X_plain.cols();
    Xj.resize(rows, cols);
    if (party_id == 2)
    {
    }
    else
    {
        for (size_t i = 0; i < rows; i++)
            for (size_t k = 0; k < cols; k++)
            {
                int64_t temp = prgs[1 - party_id].pop_int64();
                Xj(i, k) = X_plain(i, k) * party_id;
                Xj(i, k) = party_id == 0 ? Xj(i, k) - temp : Xj(i, k) + temp;
            }
    }
}

void MyEngine::reveal(const VecI64 &Xj, VecI64 &X_plain)
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
        VecI64 X0, X1;
        string2vecI64(message_recv[0], X0);
        string2vecI64(message_recv[1], X1);
        X_plain = X0 + X1;
    }
    else
    {

        vecI642string(Xj, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2vecI64(message_recv[1 - party_id], X_plain);
        X_plain = X_plain + Xj;
    }
}

void MyEngine::reveal_matrix(const MatrixXl &Xj, MatrixXl &X_plain)
{
    size_t rows = Xj.rows(), cols = Xj.cols();
    X_plain.resize(rows, cols);
    if (party_id == 2)
    {
        for (size_t i = 0; i < 2; i++)
        {
            message_send[i].clear();
        }
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        MatrixXl X0(rows, cols), X1(rows, cols);
        string2matrix(message_recv[0], X0);
        string2matrix(message_recv[1], X1);
        X_plain = X0 + X1;
    }
    else
    {

        matrix2string(Xj, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2matrix(message_recv[1 - party_id], X_plain);
        X_plain = X_plain + Xj;
    }
}

void MyEngine::matrix_mul(const MatrixXl &Xj, const MatrixXl &Yj, MatrixXl &Zj)
{
    size_t rows = Xj.rows();
    size_t lens = Xj.cols();
    assert(Yj.rows() == lens);
    size_t cols = Yj.cols();
    Zj.resize(rows, cols);
    if (party_id == 2)
    {
        MatrixXl A[2];
        MatrixXl B[2];
        MatrixXl C[2];
        for (size_t j = 0; j < 2; j++)
        {
            A[j].resize(rows, lens);
            for (size_t k = 0; k < rows; k++)
                for (size_t l = 0; l < lens; l++)
                    A[j](k, l) = prgs[j].int_buff[k * lens + l];
        }
        size_t temp_index = rows * lens;
        for (size_t j = 0; j < 2; j++)
        {
            B[j].resize(lens, cols);
            for (size_t k = 0; k < lens; k++)
                for (size_t l = 0; l < cols; l++)
                    B[j](k, l) = prgs[j].int_buff[temp_index + k * cols + l];
        }

        C[0].resize(rows, cols);
        for (size_t i = 0; i < rows; i++)
            for (size_t k = 0; k < cols; k++)
                C[0](i, k) = prgs[2].int_buff[i * cols + k];

        C[1] = (A[0] + A[1]) * (B[0] + B[1]) - C[0];
        // assert(C[1].rows() == rows && C[1].cols() == cols);
        for (size_t j = 0; j < 2; j++)
            matrix2string(C[j], message_send[j]);

        messenger->send_and_recv(message_send, message_recv);
    }
    else
    {
        MatrixXl Aj(rows, lens), Bj(lens, cols), Cj(rows, cols);
        Aj.resize(rows, lens);
        for (size_t k = 0; k < rows; k++)
            for (size_t l = 0; l < lens; l++)
                Aj(k, l) = prgs[2].int_buff[k * lens + l];

        size_t temp_index = rows * lens;

        Bj.resize(lens, cols);
        for (size_t k = 0; k < lens; k++)
            for (size_t l = 0; l < cols; l++)
                Bj(k, l) = prgs[2].int_buff[temp_index + k * cols + l];

        string temp1, temp2;
        MatrixXl Ej = Xj - Aj;
        MatrixXl Fj = Yj - Bj;
        matrix2string(Ej, temp1);
        matrix2string(Fj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear();
        messenger->send_and_recv(message_send, message_recv);
        string2matrix(message_recv[2], Cj);
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

void MyEngine::vector_mul(const VecI64 &Xj, const VecI64 &Yj, VecI64 &Zj)
{
    size_t vec_size = Xj.size();
    assert(Yj.size() == vec_size);
    Zj.resize(vec_size);
    if (party_id == 2)
    {
        VecI64 A[2];
        VecI64 B[2];
        VecI64 C[2];
        for (size_t j = 0; j < 2; j++)
        {
            A[j].resize(vec_size);
            B[j].resize(vec_size);
            C[j].resize(vec_size);
            for (size_t k = 0; k < vec_size; k++)
            {
                A[j][k] = prgs[j].int_buff[k];
                B[j][k] = prgs[j].int_buff[k + vec_size];
            }
        }
        for (size_t k = 0; k < vec_size; k++)
        {
            C[0][k] = prgs[2].int_buff[k];
            C[1][k] = (A[0][k] + A[1][k]) * (B[0][k] + B[1][k]) - C[0][k];
        }

        for (size_t j = 0; j < 2; j++)
            vecI642string(C[j], message_send[j]);

        messenger->send_and_recv(message_send, message_recv);
    }
    else
    {
        VecI64 Aj(vec_size), Bj(vec_size), Cj(vec_size), Ej(vec_size), Fj(vec_size);
        for (size_t k = 0; k < vec_size; k++)
        {
            Aj[k] = prgs[2].int_buff[k];
            Bj[k] = prgs[2].int_buff[k + vec_size];
        }
        string temp1, temp2;
        Ej = Xj - Aj;
        Fj = Yj - Bj;
        vecI642string(Ej, temp1);
        vecI642string(Fj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear();
        messenger->send_and_recv(message_send, message_recv);
        string2vecI64(message_recv[2], Cj);
        temp1 = message_recv[1 - party_id].substr(0, temp1.size());
        temp2 = message_recv[1 - party_id].substr(temp1.size(), temp2.size());
        VecI64 E_plain(vec_size);
        VecI64 F_plain(vec_size);
        string2vecI64(temp1, E_plain);
        string2vecI64(temp2, F_plain);
        E_plain = E_plain + Ej;
        F_plain = F_plain + Fj;
        Zj = Cj + (E_plain * Bj) + (Aj * F_plain) + (E_plain * F_plain * (int64_t)party_id);
    }
}

/// @brief X R 都在[0, L-2]中
/// @param Xbits_j
/// @param Beta1_plain
/// @param R_palin
/// @param Beta2_plain
void MyEngine::privateCompare(const vector<zpshare> &Xbits_j, const vecbool &Beta1_plain, const VecI64 &R_palin, vecbool &Beta2_plain)
{
    size_t vec_size = Xbits_j.size();
    assert(Beta1_plain.size() == vec_size);
    assert(R_palin.size() == vec_size);
    Beta2_plain.resize(vec_size);
    if (party_id == 0 || party_id == 1)
    {
        // P0 P1公共随机数
        vector<zpshare> S_plain;
        S_plain.resize(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            for (size_t j = 0; j < 64; j++)
            {
                S_plain[i][j] = prgs[1 - party_id].int_buff[i * 64 + j];
                S_plain[i][j] = S_plain[i][j] % P;
                if (0 == S_plain[i][j])
                {
                    S_plain[i][j] = 1;
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
                else
                {
                    Cbits_j[i][k] = (P - 1) * party_id * Tbits_plain[i][k] + Xbits_j[i][k] + (uint16_t)party_id + sum_w;
                    Cbits_j[i][k] = Cbits_j[i][k] % P;
                    Wbits_j[i][k] = Xbits_j[i][k] + (uint16_t)party_id * Tbits_plain[i][k] + (P - 2) * Tbits_plain[i][k] * Xbits_j[i][k];
                    Wbits_j[i][k] = Wbits_j[i][k] % P;
                    sum_w = sum_w + Wbits_j[i][k];
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
                else
                {
                    assert(Wbits_debug[i][k] == Xbits_debug[i][k] ^ Tbits_plain[i][k]);
                    assert(Cbits_debug[i][k] == Xbits_debug[i][k] + 1 - Tbits_plain[i][k] + sum_w);
                    sum_w = sum_w + Wbits_debug[i][k];
                }
            }
        }

#endif
    }
    else // party_id == 2
    {
        for (size_t i = 0; i < 2; i++)
        {
            message_send[i].clear();
        }
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

/// @brief 计算msb（2*A），2 * A 不会发生有符号数算术溢出时，msb（2*A）= msb（A）
/// @param Aj
/// @param Bj
void MyEngine::msb(const VecI64 &Aj, vecbool &Bj)
{
    size_t vec_size = Aj.size();
    Bj.resize(vec_size);
    VecI64 Dj = 2 * Aj;
    // 计算 msb(Dj)
    VecI64 Yj(vec_size);
    for (size_t k = 0; k < vec_size; k++)
    {
        Yj[k] = add_L_1(Dj[k], Dj[k]);
    }
    // msb(Dj) == lsb(Yj)
#if defined MYDEBUG
    VecI64 A_debug(vec_size), D_debug(vec_size), D_debug_L_1(vec_size), Y_debug(vec_size);
    reveal(Aj, A_debug);
    reveal(Dj, D_debug);
    reveal_L_1(Dj, D_debug_L_1);
    reveal_L_1(Yj, Y_debug);
    for (size_t i = 0; i < vec_size; i++)
    {
        assert((A_debug[i] < 1l << 62) && (A_debug[i] >= -1l << 62));
        assert((D_debug[i] < 0) == (A_debug[i] < 0));

        assert((D_debug[i] < 0) == (D_debug_L_1[i] < 0));

        assert((D_debug_L_1[i] < 0) == ((uint64_t)Y_debug[i] % 2));
    }
#endif
    if (party_id == 2)
    {
        // X_share = Xj
        VecI64 X_plain(vec_size), X_share[2];
        vecbool Xlsb[2];
        vector<zpshare> Xbits_plain(vec_size), X_bits[2];
        for (size_t j = 0; j < 2; j++)
        {
            X_share[j].resize(vec_size);
            Xlsb[j].resize(vec_size);
            X_bits[j].resize(vec_size);
            for (size_t k = 0; k < vec_size; k++)
            {
                X_share[j][k] = prgs[j].int_buff[k] % L_1;
            }
        }
        for (size_t k = 0; k < vec_size; k++)
        {
            X_plain[k] = add_L_1(X_share[0][k], X_share[1][k]);
            bitset<64> btst_x(X_plain[k]);
            Xlsb[0][k] = this->prgs[2].bool_buff[k];
            Xlsb[1][k] = (bool)btst_x[0] ^ Xlsb[0][k];
            for (size_t l = 0; l < 64; l++)
            {
                Xbits_plain[k][l] = (bool)btst_x[l];
                X_bits[0][k][l] = this->prgs[2].int_buff[k * 64 + l];
                X_bits[0][k][l] = X_bits[0][k][l] % P;

                X_bits[1][k][l] = Xbits_plain[k][l] + P - X_bits[0][k][l];
                X_bits[1][k][l] = X_bits[1][k][l] % P;
            }
        }
        for (size_t j = 0; j < 2; j++)
        {
            string temp1, temp2;
            vecbool2string(Xlsb[j], temp1);
            vecZpshare2string(X_bits[j], temp2);
            message_send[j] = temp1 + temp2;
        }
        messenger->send_and_recv(message_send, message_recv); // 一轮通信
        vecbool Beta1_plain(vec_size), dummy_Beta0(vec_size);
        privateCompare(Xbits_plain, dummy_Beta0, X_plain, Beta1_plain); // 二轮通信
        vecbool Beta1_share[2];
        for (size_t j = 0; j < 2; j++)
        {
            Beta1_share[j].resize(vec_size);
        }
        for (size_t k = 0; k < vec_size; k++)
        {
            Beta1_share[0][k] = prgs[2].bool_buff[vec_size + k];
            Beta1_share[1][k] = Beta1_plain[k] ^ Beta1_share[0][k];
        }
        for (size_t j = 0; j < 2; j++)
        {
            vecbool2string(Beta1_share[j], message_send[j]);
        }
        messenger->send_and_recv(message_send, message_recv); // 三轮通信
#if defined MYDEBUG
        VecI64 X_debug(vec_size);
        vecbool Beta_debug(vec_size), B_debug(vec_size), Xlsb_debug(vec_size);
        reveal_L_1(X_plain, X_debug);
        reveal_vec_bools(Beta1_plain, Beta_debug);
        reveal_vec_bools(Beta1_plain, B_debug);
        reveal_vec_bools(Xlsb[0], Xlsb_debug);

#endif
    }
    else
    {
        VecI64 Xj(vec_size);
        for (size_t k = 0; k < vec_size; k++)
            Xj[k] = prgs[2].int_buff[k] % L_1;
        vecbool Beta0_plain(vec_size);
        for (size_t k = 0; k < vec_size; k++)
            Beta0_plain[k] = prgs[1 - party_id].bool_buff[k];
        VecI64 Rj(vec_size);
        for (size_t k = 0; k < vec_size; k++)
        {
            Rj[k] = add_L_1(Xj[k], Yj[k]);
        }

        message_send[2].clear();
        vecI642string(Rj, message_send[1 - party_id]);
        messenger->send_and_recv(message_send, message_recv); // 一轮通信
        VecI64 R_plain;
        string2vecI64(message_recv[1 - party_id], R_plain);
        for (size_t k = 0; k < vec_size; k++)
        {
            R_plain[k] = add_L_1(R_plain[k], Rj[k]);
        }

        vecbool Xlsb_j(vec_size);
        vector<zpshare> Xbits_j(vec_size);
        string temp1, temp2;
        temp1 = message_recv[2].substr(0, vec_size);
        temp2 = message_recv[2].substr(vec_size, vec_size * 64 * sizeof(uint16_t));
        string2vecbool(temp1, Xlsb_j);
        string2vecZpshare(temp2, Xbits_j);
        vecbool dummy_Beta1(vec_size);
        privateCompare(Xbits_j, Beta0_plain, R_plain, dummy_Beta1); // 二轮通信
        message_send[2].clear();
        message_send[1 - party_id].clear();
        messenger->send_and_recv(message_send, message_recv); // 三轮通信
        vecbool Beta1j(vec_size);
        string2vecbool(message_recv[2], Beta1j);
        // Beta_j
        for (size_t k = 0; k < vec_size; k++)
        {
            dummy_Beta1[k] = party_id == 0 ? Beta1j[k] : Beta1j[k] ^ Beta0_plain[k];
            bool temp_R0 = (uint64_t)R_plain[k] % 2;
            Bj[k] = party_id == 1 ? dummy_Beta1[k] ^ Xlsb_j[k] : dummy_Beta1[k] ^ Xlsb_j[k] ^ temp_R0;
        }
#if defined MYDEBUG
        VecI64 X_debug(vec_size);
        vecbool Beta_debug(vec_size), B_debug(vec_size), Xlsb_debug(vec_size);
        reveal_L_1(Xj, X_debug);
        reveal_vec_bools(dummy_Beta1, Beta_debug);
        reveal_vec_bools(Bj, B_debug);
        reveal_vec_bools(Xlsb_j, Xlsb_debug);
        for (size_t i = 0; i < vec_size; i++)
        {
            assert(Beta_debug[i] == ((uint64_t)X_debug[i] > (uint64_t)R_plain[i]));
            assert(B_debug[i] == (Beta_debug[i] ^ Xlsb_debug[i] ^ ((uint64_t)R_plain[i] % 2)));
            assert(B_debug[i] == ((uint64_t)Y_debug[i] % 2));
        }
#endif
    }
}

void MyEngine::bit_inject(const VecI64 &Xj, const vecbool &Aj, VecI64 &Zj)
{
    size_t vec_size = Xj.size();
    assert(Aj.size() == vec_size);
    Zj.resize(vec_size);
    if (party_id == 1 || party_id == 0)
    {
        vecbool Bj(vec_size);
        VecI64 Uj(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Uj[i] = prgs[2].int_buff[i];
            Bj[i] = prgs[2].bool_buff[i];
        }
        VecI64 Yj(vec_size);
        vecbool Cj(vec_size);
        for (size_t i = 0; i < vec_size; i++)
        {
            Cj[i] = Bj[i] ^ Aj[i];
        }
        Yj = Xj - Uj;
        string temp1, temp2;
        vecI642string(Yj, temp1);
        vecbool2string(Cj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear();
        messenger->send_and_recv(message_send, message_recv); // 一轮通信
        VecI64 Y_plain(vec_size);
        vecbool C_plain(vec_size);
        temp1 = message_recv[1 - party_id].substr(0, vec_size * sizeof(int64_t));
        temp2 = message_recv[1 - party_id].substr(vec_size * sizeof(int64_t), vec_size);
        string2vecI64(temp1, Y_plain);
        string2vecbool(temp2, C_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            Y_plain[i] = Y_plain[i] + Yj[i];
            C_plain[i] = C_plain[i] ^ Cj[i];
        }
        VecI64 Vj(vec_size), Wj(vec_size), Sj(vec_size), Tj(vec_size);
        {
            string str_V, str_W, str_S, str_T;
            size_t bat_size = vec_size * sizeof(int64_t);
            str_V = message_recv[2].substr(0, bat_size);
            str_W = message_recv[2].substr(1 * bat_size, bat_size);
            str_S = message_recv[2].substr(2 * bat_size, bat_size);
            str_T = message_recv[2].substr(3 * bat_size, bat_size);
            string2vecI64(str_V, Vj);
            string2vecI64(str_W, Wj);
            string2vecI64(str_S, Sj);
            string2vecI64(str_T, Tj);
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            if (C_plain[i] == 0) // B[i] == A[i]
            {
                Zj[i] = Vj[i] + Y_plain[i] * (Bj[i] - 2 * Sj[i]);
            }
            else // A[i] == B[i] ^ 1
            {
                Zj[i] = Wj[i] + Y_plain[i] * ((party_id ^ Bj[i]) - 2 * Tj[i]);
            }
        }
#if defined MYDEBUG
        vecbool B_debug(vec_size), A_debug(vec_size);
        reveal_vec_bools(Bj, B_debug);
        reveal_vec_bools(Aj, A_debug);
        VecI64 X_debug(vec_size), U_debug(vec_size);
        reveal(Xj, X_debug);
        reveal(Uj, U_debug);
        VecI64 V_debug(vec_size), W_debug(vec_size),S_debug(vec_size), T_debug(vec_size);
        reveal(Vj, V_debug);
        reveal(Wj, W_debug);
        reveal(Sj, S_debug);
        reveal(Tj, T_debug);

        // printf("Bj\n");
        // printVecbool(Bj);
        // printf("B_debug\n");
        // printVecbool(B_debug);
        // printf("A_debug\n");
        // printVecbool(A_debug);
        // printf("C_plain\n");
        // printVecbool(C_plain);

        // printf("X_debug\n");
        // printVecI64(X_debug);
        // printf("U_debug\n");
        // printVecI64(U_debug);
        // printf("Y_plain\n");
        // printVecI64(Y_plain);

        // printf("V_debug\n");
        // printVecI64(V_debug);
        // printf("W_debug\n");
        // printVecI64(W_debug);
        // printf("S_debug\n");
        // printVecI64(S_debug);
        // printf("T_debug\n");
        // printVecI64(T_debug);

#endif
    }
    else
    {
        vecbool B[2], B_plain(vec_size);
        VecI64 U[2], U_plain(vec_size);
        for (size_t j = 0; j < 2; j++)
        {
            B[j].resize(vec_size);
            U[j].resize(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                B[j][i] = prgs[j].bool_buff[i];
                U[j][i] = prgs[j].int_buff[i];
            }
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            B_plain[i] = B[0][i] ^ B[1][i];
            U_plain[i] = U[0][i] + U[1][i];
        }
        VecI64 V[2], W[2], S[2], T[2];
        for (size_t j = 0; j < 2; j++)
        {
            V[j].resize(vec_size);
            W[j].resize(vec_size);
            S[j].resize(vec_size);
            T[j].resize(vec_size);
        }
        for (size_t i = 0; i < vec_size; i++)
        {
            V[0][i] = prgs[2].int_buff[i];
            V[1][i] = B_plain[i] * U_plain[i] - V[0][i];

            W[0][i] = prgs[2].int_buff[vec_size + i];
            W[1][i] = (1 ^ B_plain[i]) * U_plain[i] - W[0][i];

            S[0][i] = prgs[2].int_buff[2 * vec_size + i];
            S[1][i] = B[0][i] * B[1][i] - S[0][i];

            T[0][i] = prgs[2].int_buff[3 * vec_size + i];
            T[1][i] = B[0][i] * (1 ^ B[1][i]) - T[0][i];
        }
        for (size_t j = 0; j < 2; j++)
        {
            string str_V, str_W, str_S, str_T;
            vecI642string(V[j], str_V);
            vecI642string(W[j], str_W);
            vecI642string(S[j], str_S);
            vecI642string(T[j], str_T);
            message_send[j] = str_V + str_W + str_S + str_T;
        }
        messenger->send_and_recv(message_send, message_recv); // 一轮通信
#if defined MYDEBUG
        vecbool B_debug(vec_size), A_debug(vec_size);
        reveal_vec_bools(B_plain, B_debug);
        reveal_vec_bools(Aj, A_debug);
        VecI64 X_debug(vec_size), U_debug(vec_size);
        reveal(Xj, X_debug);
        reveal(U_plain, U_debug);
        VecI64 V_debug(vec_size), W_debug(vec_size),S_debug(vec_size), T_debug(vec_size);
        reveal(V[0], V_debug);
        reveal(W[0], W_debug);
        reveal(S[0], S_debug);
        reveal(T[0], T_debug);

        // printf("B_plain\n");
        // printVecbool(B_plain);

        // printf("U_plain\n");
        // printVecI64(U_plain);
#endif
    }
}

/// @brief 和securenn中的不一样，Xj分别被P0 P1持有
/// @param Xj
/// @param X_plain
void MyEngine::reveal_vec_bools(const vecbool &Xj, vecbool &X_plain)
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
        vecbool X0, X1;
        string2vecbool(message_recv[0], X0);
        string2vecbool(message_recv[1], X1);
        for (size_t k = 0; k < vec_size; k++)
        {
            X_plain[k] = X0[k] ^ X1[k];
        }
    }
    else
    {
        vecbool2string(Xj, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2vecbool(message_recv[1 - party_id], X_plain);
        for (size_t k = 0; k < vec_size; k++)
        {
            X_plain[k] = X_plain[k] ^ Xj[k];
        }
    }
}

void MyEngine::reveal_L_1(const VecI64 &Xj, VecI64 &X_plain)
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
        VecI64 X0, X1;
        string2vecI64(message_recv[0], X0);
        string2vecI64(message_recv[1], X1);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = add_L_1(X0[i], X1[i]);
        }
    }
    else
    {
        vecI642string(Xj, message_send[1 - party_id]);
        message_send[2] = message_send[1 - party_id];
        messenger->send_and_recv(message_send, message_recv); // 第一轮通信
        string2vecI64(message_recv[1 - party_id], X_plain);
        for (size_t i = 0; i < vec_size; i++)
        {
            X_plain[i] = add_L_1(Xj[i], X_plain[i]);
        }
    }
}

void MyEngine::reveal_zp(const vector<zpshare> &Xbits_j, vector<zpshare> &Xbits_plain)
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
