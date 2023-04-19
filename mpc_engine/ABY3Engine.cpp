#include "ABY3Engine.hpp"

// #define MYDEBUG

static VecI64 bitwise_and(const VecI64 &a, const VecI64 &b)
{
    size_t vec_size = a.size();
    assert(b.size() == vec_size);
    VecI64 ret(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        ret[i] = a[i] & b[i];
    }
    return ret;
}

static VecI64 bitwise_xor(const VecI64 &a, const VecI64 &b)
{
    size_t vec_size = a.size();
    assert(b.size() == vec_size);
    VecI64 ret(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        ret[i] = a[i] ^ b[i];
    }
    return ret;
}

ABY3Engine::ABY3Engine(VirtuleMessenger *_messenger, size_t _party_id) : messenger(_messenger), party_id(_party_id)
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

ABY3Engine::~ABY3Engine()
{
}

/// @brief
/// @param X_plain Xj[0] = X_j Xj[1] = X_{nextParty}
/// @param Xj
void ABY3Engine::share(const VecI64 &X_plain, array<VecI64, 2> &Xj)
{
    size_t vec_size = X_plain.size();
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;

    Xj[0].resize(vec_size);
    Xj[1].resize(vec_size);
    VecI64 Rj(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        Rj(i) = prgs[nextParty].pop_int64() - prgs[prevParty].pop_int64();
    }
    Xj[0] = party_id == 0 ? Rj + X_plain : Rj;
    message_send[nextParty].clear();
    vecI642string(Xj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecI64(message_recv[nextParty], Xj[1]);
}

void ABY3Engine::share_matrix(const MatrixXl &X_plain, array<MatrixXl, 2> &Xj)
{
    size_t rows = X_plain.rows();
    size_t cols = X_plain.cols();

    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    Xj[0].resize(rows, cols);
    Xj[1].resize(rows, cols);

    MatrixXl Rj(rows, cols);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            Rj(i, j) = prgs[nextParty].pop_int64() - prgs[prevParty].pop_int64();
        }
    }
    Xj[0] = party_id == 0 ? Rj + X_plain : Rj;
    message_send[nextParty].clear();
    matrix2string(Xj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2matrix(message_recv[nextParty], Xj[1]);
}

void ABY3Engine::share_bools(const vecbool &X_plain, array<vecbool, 2> &Xj)
{
    size_t vec_size = X_plain.size();
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;

    Xj[0].resize(vec_size);
    Xj[1].resize(vec_size);
    vecbool Rj(vec_size);
    for (size_t i = 0; i < vec_size; i++)
    {
        Rj[i] = prgs[nextParty].pop_bool() ^ prgs[prevParty].pop_bool();
        Xj[0][i] = party_id == 0 ? Rj[i] ^ X_plain[i] : Rj[i];
    }

    message_send[nextParty].clear();
    vecbool2string(Xj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecbool(message_recv[nextParty], Xj[1]);
}

/// @brief
/// @param X_plain Xj[0] = X_j Xj[1] = X_{nextParty}
/// @param Xj
void ABY3Engine::reveal(const array<VecI64, 2> &Xj, VecI64 &X_plain)
{
    size_t vec_size = X_plain.size();
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    message_send[prevParty].clear();
    vecI642string(Xj[0], message_send[nextParty]);
    messenger->send_and_recv(message_send, message_recv);
    VecI64 X_prev(vec_size);
    string2vecI64(message_recv[prevParty], X_prev);
    X_plain = X_prev + Xj[0] + Xj[1];
}

void ABY3Engine::reveal_matrix(const array<MatrixXl, 2> &Xj, MatrixXl &X_plain)
{
    size_t rows = Xj[0].rows();
    size_t cols = Xj[0].cols();
    assert(Xj[1].rows() == rows && Xj[1].cols() == cols);
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    message_send[prevParty].clear();
    matrix2string(Xj[0], message_send[nextParty]);
    messenger->send_and_recv(message_send, message_recv);
    MatrixXl X_prev(rows, cols);
    string2matrix(message_recv[prevParty], X_prev);
    X_plain = X_prev + Xj[0] + Xj[1];
}

void ABY3Engine::reveal_bools(const array<vecbool, 2> &Xj, vecbool &X_plain)
{
    size_t vec_size = Xj[0].size();
    X_plain.resize(vec_size);
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    message_send[prevParty].clear();
    vecbool2string(Xj[0], message_send[nextParty]);
    messenger->send_and_recv(message_send, message_recv);
    vecbool X_prev(vec_size);
    string2vecbool(message_recv[prevParty], X_prev);
    for (size_t i = 0; i < vec_size; i++)
    {
        X_plain[i] = Xj[0][i] ^ Xj[1][i] ^ X_prev[i];
    }
}

void ABY3Engine::matrix_mul(const array<MatrixXl, 2> &Xj, const array<MatrixXl, 2> &Yj, array<MatrixXl, 2> &Zj)
{
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    size_t rows = Xj[0].rows();
    size_t lens = Xj[0].cols();
    assert(Yj[0].rows() == lens);
    size_t cols = Yj[0].cols();
    assert(Xj[1].rows() == rows && Xj[1].cols() == lens);
    assert(Yj[1].rows() == lens && Yj[1].cols() == cols);
    Zj[0].resize(rows, cols);
    Zj[1].resize(rows, cols);
    Zj[0] = Xj[0] * Yj[0] + Xj[0] * Yj[1] + Xj[1] * Yj[0];
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t k = 0; k < cols; k++)
        {
            int64_t temp = prgs[nextParty].pop_int64() - prgs[prevParty].pop_int64();
            Zj[0](i, k) = Zj[0](i, k) + temp;
        }
    }
    message_send[nextParty].clear();
    matrix2string(Zj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2matrix(message_recv[nextParty], Zj[1]);
#if defined MYDEBUG
    MatrixXl X_debug(rows, lens), Y_debug(lens, cols), Z_debug(rows, cols);
    reveal_matrix(Xj, X_debug);
    reveal_matrix(Yj, Y_debug);
    reveal_matrix(Zj, Z_debug);

    assert(Z_debug == X_debug * Y_debug);
#endif
}

void ABY3Engine::vector_mul(const array<VecI64, 2> &Xj, const array<VecI64, 2> &Yj, array<VecI64, 2> &Zj)
{
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    size_t vec_size = Xj[0].size();
    assert(Xj[1].size() == vec_size && Yj[0].size() == vec_size && Yj[1].size() == vec_size);
    Zj[0].resize(vec_size);
    Zj[1].resize(vec_size);
    Zj[0] = Xj[0] * Yj[0] + Xj[0] * Yj[1] + Xj[1] * Yj[0];
    for (size_t i = 0; i < vec_size; i++)
    {
        int64_t temp = prgs[nextParty].pop_int64() - prgs[prevParty].pop_int64();
        Zj[0](i) = Zj[0](i) + temp;
    }
    message_send[nextParty].clear();
    vecI642string(Zj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecI64(message_recv[nextParty], Zj[1]);

#if defined MYDEBUG
    VecI64 X_debug, Y_debug, Z_debug;
    reveal(Xj, X_debug);
    reveal(Yj, Y_debug);
    reveal(Zj, Z_debug);
    for (size_t i = 0; i < vec_size; i++)
    {
        assert(Z_debug[i] == X_debug[i] * Y_debug[i]);
    }

#endif
}

void ABY3Engine::vecbool_mul(const array<vecbool, 2> &Xj, const array<vecbool, 2> &Yj, array<vecbool, 2> &Zj)
{
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    size_t vec_size = Xj[0].size();
    assert(Xj[1].size() == vec_size && Yj[0].size() == vec_size && Yj[1].size() == vec_size);
    Zj[0].resize(vec_size);
    Zj[1].resize(vec_size);

    for (size_t i = 0; i < vec_size; i++)
    {
        Zj[0][i] = (Xj[0][i] & Yj[0][i]) ^ (Xj[1][i] & Yj[0][i]) ^ (Xj[0][i] & Yj[1][i]);
    }
    message_send[nextParty].clear();
    vecbool2string(Zj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecbool(message_recv[nextParty], Zj[1]);
#if defined MYDEBUG
    vecbool X_debug, Y_debug, Z_debug;
    reveal_bools(Xj, X_debug);
    reveal_bools(Yj, Y_debug);
    reveal_bools(Zj, Z_debug);
    for (size_t i = 0; i < vec_size; i++)
    {
        assert(Z_debug[i] == X_debug[i] * Y_debug[i]);
    }

#endif
}

/// @brief 比特秘密分享的暂时形式
/// @param Xj
/// @param X_plain
void ABY3Engine::reveal_bits(const array<VecI64, 2> &Xj, VecI64 &X_plain)
{
    size_t vec_size = X_plain.size();
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    message_send[prevParty].clear();
    vecI642string(Xj[0], message_send[nextParty]);
    messenger->send_and_recv(message_send, message_recv);
    VecI64 X_prev(vec_size);
    string2vecI64(message_recv[prevParty], X_prev);
    VecI64 temp = bitwise_xor(Xj[0], Xj[1]);
    X_plain = bitwise_xor(X_prev, temp);
}

/// @brief P = _2C ^ S; G = _2C & S
/// @param Xj
/// @param Gj
/// @param Pj
void ABY3Engine::add_prepare(const array<VecI64, 2> &Xj, array<VecI64, 2> &Gj, array<VecI64, 2> &Pj)
{
    size_t vec_size = Xj[0].size();
    assert(Xj[1].size() == vec_size);
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    // 比特秘密分享 S = S0 ^ S1 ^ S2
    array<VecI64, 2> Cj, Sj;
    Sj[0] = Xj[0];
    Sj[1] = Xj[1];
    Cj[0] = bitwise_and(Xj[0], Xj[1]);
    message_send[nextParty].clear();
    vecI642string(Cj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecI64(message_recv[nextParty], Cj[1]);
    array<VecI64, 2> _2Cj;
    _2Cj[0] = Cj[0] * 2;
    _2Cj[1] = Cj[1] * 2;

    Pj[0] = bitwise_xor(_2Cj[0], Sj[0]);
    Pj[1] = bitwise_xor(_2Cj[1], Sj[1]);
    VecI64 temp1, temp2, temp3;
    temp1 = bitwise_and(_2Cj[0], Sj[0]);
    temp2 = bitwise_and(_2Cj[0], Sj[1]);
    temp3 = bitwise_and(_2Cj[1], Sj[0]);

    Gj[0] = bitwise_xor(temp1, temp2);
    Gj[0] = bitwise_xor(Gj[0], temp3);
    message_send[nextParty].clear();
    vecI642string(Gj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecI64(message_recv[nextParty], Gj[1]);
#if defined MYDEBUG
    VecI64 X_debug, _2C_debug, S_debug;
    reveal(Xj, X_debug);
    reveal_bits(_2Cj, _2C_debug);
    reveal_bits(Sj, S_debug);

    for (size_t i = 0; i < vec_size; i++)
    {
        assert(S_debug[i] + _2C_debug[i] == X_debug[i]);
    }
    VecI64 G_debug, P_debug;
    reveal_bits(Gj, G_debug);
    reveal_bits(Pj, P_debug);
    VecI64 temp = bitwise_and(S_debug, _2C_debug);

    for (size_t i = 0; i < vec_size; i++)
    {
        assert(P_debug[i] == S_debug[i] ^ _2C_debug[i]);
        assert(G_debug[i] == temp[i]);
        assert(temp[i] == (S_debug[i] & _2C_debug[i]));
    }
    // printf("X_debug\n");
    // printVecI64(X_debug);
    // printf("_2C_debug\n");
    // printVecI64(_2C_debug);
    // printf("S_debug\n");
    // printVecI64(S_debug);
    // printf("P_debug\n");
    // printVecI64(P_debug);
    // printf("G_debug\n");
    // printVecI64(G_debug);
#endif
}

void ABY3Engine::msb(const array<VecI64, 2> &Gj, const array<VecI64, 2> &Pj, array<vecbool, 2> &Zj)
{
    size_t vec_size = Gj[0].size();
    assert(vec_size == Gj[1].size());
    assert(vec_size == Pj[0].size());
    assert(vec_size == Pj[1].size());

    array<vecbool, 2> Gj_vec, Pj_vec;
    array<vecbool, 2> backup;
    for (size_t i = 0; i < 2; i++)
    {
        Gj_vec[i].resize(vec_size * 64);
        Pj_vec[i].resize(vec_size * 64);
        backup[i].resize(vec_size);
        for (size_t k = 0; k < vec_size; k++)
        {

            bitset<64> tempG(Gj[i][k]);
            bitset<64> tempP(Pj[i][k]);
            for (size_t l = 0; l < 64; l++)
            {
                Gj_vec[i][k * 64 + l] = (bool)tempG[l];
                Pj_vec[i][k * 64 + l] = (bool)tempP[l];
            }
            backup[i][k] = Pj_vec[i][k * 64 + 63];
            Gj_vec[i][k * 64 + 63] = 0;
            Pj_vec[i][k * 64 + 63] = 1;
        }
    }

    array<vecbool, 2> Gj_vec_new, Pj_vec_new;
    array<vecbool, 2> op1, op2, result;

    // 6轮 temp = 32, 16, 8, 4, 2, 1;
    for (size_t temp = 32; 0 != temp; temp = temp / 2)
    {
        for (size_t i = 0; i < 2; i++)
        {
            op1[i].resize(vec_size * 2 * temp);
            op2[i].resize(vec_size * 2 * temp);
            result[i].resize(vec_size * 2 * temp);
            for (size_t k = 0; k < vec_size; k++)
            {
                for (size_t m = 0; m < temp; m++)
                {
                    op1[i][k * 2 * temp + 2 * m] = Pj_vec[i][k * 2 * temp + 2 * m + 1];
                    op2[i][k * 2 * temp + 2 * m] = Gj_vec[i][k * 2 * temp + 2 * m];

                    op1[i][k * 2 * temp + 2 * m + 1] = Pj_vec[i][k * 2 * temp + 2 * m + 1];
                    op2[i][k * 2 * temp + 2 * m + 1] = Pj_vec[i][k * 2 * temp + 2 * m];
                }
            }
        }
        vecbool_mul(op1, op2, result);

        for (size_t i = 0; i < 2; i++)
        {
            Gj_vec_new[i].resize(vec_size * temp);
            Pj_vec_new[i].resize(vec_size * temp);
            for (size_t k = 0; k < vec_size; k++)
            {
                for (size_t m = 0; m < temp; m++)
                {
                    Gj_vec_new[i][k * temp + m] = result[i][k * 2 * temp + 2 * m] ^ Gj_vec[i][k * 2 * temp + 2 * m + 1];
                    Pj_vec_new[i][k * temp + m] = result[i][k * 2 * temp + 2 * m + 1];
                }
            }
        }

#if defined MYDEBUG
        vecbool Gj_vec_debug, Pj_vec_debug;
        reveal_bools(Gj_vec, Gj_vec_debug);
        reveal_bools(Pj_vec, Pj_vec_debug);
        vecbool Gj_vec_new_debug, Pj_vec_new_debug;
        reveal_bools(Gj_vec_new, Gj_vec_new_debug);
        reveal_bools(Pj_vec_new, Pj_vec_new_debug);
        for (size_t k = 0; k < vec_size; k++)
        {
            for (size_t m = 0; m < temp; m++)
            {

                assert(Pj_vec_new_debug[k * temp + m] == (Pj_vec_debug[k * 2 * temp + 2 * m] & Pj_vec_debug[k * 2 * temp + 2 * m + 1]));

                int8_t temp_bool = Gj_vec_debug[k * 2 * temp + 2 * m] & Pj_vec_debug[k * 2 * temp + 2 * m + 1];
                temp_bool = temp_bool ^ Gj_vec_debug[k * 2 * temp + 2 * m + 1];
                assert(Gj_vec_new_debug[k * temp + m] == temp_bool);
            }
        }
        // printf("in msb round %lu\n", temp);
        // printf("Gj_vec_debug\n");
        // printVecbool(Gj_vec_debug);
        // printf("Pj_vec_debug\n");
        // printVecbool(Pj_vec_debug);
#endif
        Gj_vec = Gj_vec_new;
        Pj_vec = Pj_vec_new;
    }

    for (size_t i = 0; i < 2; i++)
    {
        Zj[i].resize(vec_size);
        for (size_t k = 0; k < vec_size; k++)
        {
            Zj[i][k] = Gj_vec[i][k] ^ backup[i][k];
        }
    }
}

void ABY3Engine::trunc(const array<VecI64, 2> &Xj, array<VecI64, 2> &Yj, size_t number)
{
    assert(number < 64);
    size_t vec_size = Xj[0].size();
    assert(vec_size == Xj[1].size());

    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    Yj[0].resize(vec_size);
    Yj[1].resize(vec_size);

    if (party_id == 0)
    {
        for (size_t i = 0; i < vec_size; i++)
        {
            Yj[0][i] = Xj[0][i] >> number;
        }
    }
    else if (party_id == 1)
    {

        for (size_t i = 0; i < vec_size; i++)
        {
            Yj[0][i] = Xj[0][i] + Xj[1][i];
            Yj[0][i] = Yj[0][i] >> number;
            Yj[0][i] = Yj[0][i] - prgs[2].pop_int64();
        }
    }
    else
    {
        for (size_t i = 0; i < vec_size; i++)
        {
            Yj[0][i] = prgs[1].pop_int64();
        }
    }
    message_send[nextParty].clear();
    vecI642string(Yj[0], message_send[prevParty]);
    messenger->send_and_recv(message_send, message_recv);
    string2vecI64(message_recv[nextParty], Yj[1]);
}

void ABY3Engine::bit_inject(const array<VecI64, 2> &Aj, const array<vecbool, 2> &Bj, array<VecI64, 2> &Cj)
{
#if defined MYDEBUG
    // printf("Aj[0]\n");
    // printVecI64(Aj[0]);
    // printf("Aj[1]\n");
    // printVecI64(Aj[1]);
    // printf("Bj[0]\n");
    // printVecbool(Bj[0]);
    // printf("Bj[1]\n");
    // printVecbool(Bj[1]);
#endif
    size_t vec_size = Aj[0].size();
    assert(vec_size == Aj[1].size());
    assert(vec_size == Bj[1].size());
    assert(vec_size == Bj[0].size());
    size_t nextParty = (party_id + 1) % 3;
    size_t prevParty = (party_id + 2) % 3;
    message_send[nextParty].clear();
    message_send[prevParty].clear();
    Cj[0].resize(vec_size);
    Cj[1].resize(vec_size);

    if (party_id == 0)
    {
        string temp0, temp1;
        VecI64 C01(vec_size), C02(vec_size);
        // P0发送者 通信前
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                C01[i] = prgs[1].pop_int64();
                C02[i] = prgs[2].pop_int64();
            }
            VecI64 M0_0(vec_size), M0_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                M0_0[i] = (0 ^ Bj[0][i] ^ Bj[1][i]) * Aj[0][i] - C01[i] - C02[i];
                M0_1[i] = (1 ^ Bj[0][i] ^ Bj[1][i]) * Aj[0][i] - C01[i] - C02[i];
            }
            VecI64 W02_0(vec_size), W02_1(vec_size);
            VecI64 W01_0(vec_size), W01_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W02_0[i] = prgs[2].pop_int64();
                W02_1[i] = prgs[2].pop_int64();
                W01_0[i] = prgs[1].pop_int64();
                W01_1[i] = prgs[1].pop_int64();
            }
            vecI642string(bitwise_xor(M0_0, W02_0), temp0);
            vecI642string(bitwise_xor(M0_1, W02_1), temp1);
            message_send[1] = temp0 + temp1;
            vecI642string(bitwise_xor(M0_0, W01_0), temp0);
            vecI642string(bitwise_xor(M0_1, W01_1), temp1);
            message_send[2] = temp0 + temp1;
#if defined MYDEBUG
            // printf("[LOG]:\tP0发送者 通信前\n");
            // print_plus();
            // printf("C01\n");
            // printVecI64(C01);
            // printf("C02\n");
            // printVecI64(C02);
            // printf("M0_0\n");
            // printVecI64(M0_0);
            // printf("M0_1\n");
            // printVecI64(M0_1);
            // printf("W02_0\n");
            // printVecI64(W02_0);
            // printf("W02_1\n");
            // printVecI64(W02_1);
            // printf("W01_0\n");
            // printVecI64(W01_0);
            // printf("W01_1\n");
            // printVecI64(W01_1);
            // print_minus();
#endif
        }
        // P1发送者 通信前
        VecI64 C10(vec_size);
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                C10[i] = prgs[1].pop_int64();
            }
            VecI64 W10_0(vec_size), W10_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W10_0[i] = prgs[1].pop_int64();
                W10_1[i] = prgs[1].pop_int64();
            }
            VecI64 W10_01(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W10_01[i] = Bj[0][i] == 1 ? W10_1[i] : W10_0[i];
            }
            vecI642string(W10_01, temp0);
            message_send[2] = message_send[2] + temp0;
#if defined MYDEBUG
            // printf("[LOG]:\tP1发送者 通信前\n");
            // print_plus();
            // printf("C10\n");
            // printVecI64(C10);
            // printf("W10_0\n");
            // printVecI64(W10_0);
            // printf("W10_1\n");
            // printVecI64(W10_1);
            // printf("W10_01\n");
            // printVecI64(W10_01);
            // print_minus();
#endif
        }
        messenger->send_and_recv(message_send, message_recv);
        array<VecI64, 2> A0j, A12j;
        for (size_t i = 0; i < 2; i++)
        {
            A0j[i].resize(vec_size);
            A12j[i].resize(vec_size);
        }

        // P0发送者 通信后
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                A0j[0][i] = C02[i];
                A0j[1][i] = C01[i];
            }
#if defined MYDEBUG
            // printf("[LOG]:\tP0发送者 通信后\n");
            // print_plus();
            // printf("A0j[0]\n");
            // printVecI64(A0j[0]);
            // printf("A0j[1]\n");
            // printVecI64(A0j[1]);
            // print_minus();
#endif
        }
        // P1发送者 通信后
        {

            VecI64 MW12_0(vec_size), MW12_1(vec_size);
            assert(message_recv[1].size() == 2 * vec_size * sizeof(int64_t));
            temp0 = message_recv[1].substr(0, vec_size * sizeof(int64_t));
            temp1 = message_recv[1].substr(vec_size * sizeof(int64_t), vec_size * sizeof(int64_t));
            string2vecI64(temp0, MW12_0);
            string2vecI64(temp1, MW12_1);
            VecI64 W12_01(vec_size);
            string2vecI64(message_recv[2], W12_01);
            for (size_t i = 0; i < vec_size; i++)
            {
                A12j[1][i] = C10[i];
                A12j[0][i] = Bj[0][i] == 1 ? MW12_1[i] ^ W12_01[i] : MW12_0[i] ^ W12_01[i];
            }
#if defined MYDEBUG
            // printf("[LOG]:\tP1发送者 通信后\n");
            // print_plus();

            // printf("MW12_0\n");
            // printVecI64(MW12_0);
            // printf("MW12_1\n");
            // printVecI64(MW12_1);
            // printf("W12_01\n");
            // printVecI64(W12_01);
            // printf("A12j[0]\n");
            // printVecI64(A12j[0]);
            // printf("A12j[1]\n");
            // printVecI64(A12j[1]);
            // print_minus();
#endif
        }
        for (size_t i = 0; i < 2; i++)
        {
            Cj[i] = A0j[i] + A12j[i];
        }
#if defined MYDEBUG
        vecbool B_debug;
        reveal_bools(Bj, B_debug);
        VecI64 A0_debug;
        reveal(A0j, A0_debug);
        if (party_id == 0)
            for (size_t i = 0; i < vec_size; i++)
            {
                assert(A0_debug[i] == Aj[0][i] * B_debug[i]);
            }
        VecI64 A12_debug;
        reveal(A12j, A12_debug);
        if (party_id == 1)
            for (size_t i = 0; i < vec_size; i++)
            {
                assert(A12_debug[i] == (Aj[0][i] + Aj[1][i]) * B_debug[i]);
            }
#endif
    }
    else if (party_id == 1)
    {
        string temp0, temp1;
        VecI64 C01(vec_size);
        // P0发送者 通信前
        {

            for (size_t i = 0; i < vec_size; i++)
            {
                C01[i] = prgs[0].pop_int64();
            }
            VecI64 W01_0(vec_size), W01_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W01_0[i] = prgs[0].pop_int64();
                W01_1[i] = prgs[0].pop_int64();
            }
            VecI64 W01_01(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W01_01[i] = Bj[1][i] == 1 ? W01_1[i] : W01_0[i];
            }
            vecI642string(W01_01, temp0);
            message_send[2] = temp0;
#if defined MYDEBUG
            // printf("[LOG]:\tP0发送者 通信前\n");
            // print_plus();
            // printf("C01\n");
            // printVecI64(C01);
            // printf("W01_0\n");
            // printVecI64(W01_0);
            // printf("W01_1\n");
            // printVecI64(W01_1);
            // printf("W01_01\n");
            // printVecI64(W01_01);
            // print_minus();
#endif
        }
        VecI64 C10(vec_size), C12(vec_size);
        // P1发送者 通信前
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                C10[i] = prgs[0].pop_int64();
                C12[i] = prgs[2].pop_int64();
            }
            VecI64 M1_0(vec_size), M1_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                M1_0[i] = (0 ^ Bj[0][i] ^ Bj[1][i]) * (Aj[0][i] + Aj[1][i]) - C10[i] - C12[i];
                M1_1[i] = (1 ^ Bj[0][i] ^ Bj[1][i]) * (Aj[0][i] + Aj[1][i]) - C10[i] - C12[i];
            }
            VecI64 W10_0(vec_size), W10_1(vec_size);
            VecI64 W12_0(vec_size), W12_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W10_0[i] = prgs[0].pop_int64();
                W10_1[i] = prgs[0].pop_int64();
                W12_0[i] = prgs[2].pop_int64();
                W12_1[i] = prgs[2].pop_int64();
            }
            vecI642string(bitwise_xor(M1_0, W10_0), temp0);
            vecI642string(bitwise_xor(M1_1, W10_1), temp1);
            message_send[2] = message_send[2] + temp0 + temp1;
            vecI642string(bitwise_xor(M1_0, W12_0), temp0);
            vecI642string(bitwise_xor(M1_1, W12_1), temp1);
            message_send[0] = message_send[0] + temp0 + temp1;
#if defined MYDEBUG
            // printf("[LOG]:\tP1发送者 通信前\n");
            // print_plus();
            // printf("C10\n");
            // printVecI64(C10);
            // printf("C12\n");
            // printVecI64(C12);
            // printf("M1_0\n");
            // printVecI64(M1_0);
            // printf("M1_1\n");
            // printVecI64(M1_1);
            // printf("W10_0\n");
            // printVecI64(W10_0);
            // printf("W10_1\n");
            // printVecI64(W10_1);
            // printf("W12_0\n");
            // printVecI64(W12_0);
            // printf("W12_1\n");
            // printVecI64(W12_1);
            // print_minus();
#endif
        }
        messenger->send_and_recv(message_send, message_recv);
        array<VecI64, 2> A0j, A12j;
        for (size_t i = 0; i < 2; i++)
        {
            A0j[i].resize(vec_size);
            A12j[i].resize(vec_size);
        }
        // P0发送者 通信后
        {
            VecI64 MW02_0(vec_size), MW02_1(vec_size);
            assert(message_recv[0].size() == 2 * vec_size * sizeof(int64_t));
            temp0 = message_recv[0].substr(0, vec_size * sizeof(int64_t));
            temp1 = message_recv[0].substr(vec_size * sizeof(int64_t), 2 * vec_size * sizeof(int64_t));
            string2vecI64(temp0, MW02_0);
            string2vecI64(temp1, MW02_1);
            VecI64 W02_01(vec_size);
            string2vecI64(message_recv[2], W02_01);
            for (size_t i = 0; i < vec_size; i++)
            {
                A0j[0][i] = C01[i];
                A0j[1][i] = Bj[1][i] == 1 ? MW02_1[i] ^ W02_01[i] : MW02_0[i] ^ W02_01[i];
            }
#if defined MYDEBUG
            // printf("[LOG]:\tP0发送者 通信后\n");
            // print_plus();
            // printf("MW02_0\n");
            // printVecI64(MW02_0);
            // printf("MW02_1\n");
            // printVecI64(MW02_1);
            // printf("W02_01\n");
            // printVecI64(W02_01);
            // printf("A0j[0]\n");
            // printVecI64(A0j[0]);
            // printf("A0j[1]\n");
            // printVecI64(A0j[1]);
            // print_minus();
#endif
        }
        // P1发送者 通信后
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                A12j[0][i] = C10[i];
                A12j[1][i] = C12[i];
            }
#if defined MYDEBUG
            // printf("[LOG]:\tP1发送者 通信后\n");
            // print_plus();
            // printf("A12j[0]\n");
            // printVecI64(A12j[0]);
            // printf("A12j[1]\n");
            // printVecI64(A12j[1]);
            // print_minus();
#endif
        }
        for (size_t i = 0; i < 2; i++)
        {
            Cj[i] = A0j[i] + A12j[i];
        }
#if defined MYDEBUG
        vecbool B_debug;
        reveal_bools(Bj, B_debug);
        VecI64 A0_debug;
        reveal(A0j, A0_debug);
        if (party_id == 0)
            for (size_t i = 0; i < vec_size; i++)
            {
                assert(A0_debug[i] == Aj[0][i] * B_debug[i]);
            }
        VecI64 A12_debug;
        reveal(A12j, A12_debug);
        if (party_id == 1)
            for (size_t i = 0; i < vec_size; i++)
            {
                assert(A12_debug[i] == (Aj[0][i] + Aj[1][i]) * B_debug[i]);
            }
#endif
    }
    else // P2
    {
        string temp0, temp1;
        VecI64 C02(vec_size);
        // P0发送者 通信前
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                C02[i] = prgs[0].pop_int64();
            }
            VecI64 W02_0(vec_size), W02_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W02_0[i] = prgs[0].pop_int64();
                W02_1[i] = prgs[0].pop_int64();
            }
            VecI64 W02_01(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W02_01[i] = Bj[0][i] == 1 ? W02_1[i] : W02_0[i];
            }
            vecI642string(W02_01, temp0);
            message_send[1] = temp0;
#if defined MYDEBUG
            // printf("[LOG]:\tP0发送者 通信前\n");
            // print_plus();
            // printf("C02\n");
            // printVecI64(C02);
            // printf("W02_0\n");
            // printVecI64(W02_0);
            // printf("W02_1\n");
            // printVecI64(W02_1);
            // printf("W02_01\n");
            // printVecI64(W02_01);
            // print_minus();
#endif
        }
        VecI64 C12(vec_size);
        // P1发送者 通信前
        {
            for (size_t i = 0; i < vec_size; i++)
            {
                C12[i] = prgs[1].pop_int64();
            }
            VecI64 W12_0(vec_size), W12_1(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W12_0[i] = prgs[1].pop_int64();
                W12_1[i] = prgs[1].pop_int64();
            }
            VecI64 W12_01(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                W12_01[i] = Bj[1][i] == 1 ? W12_1[i] : W12_0[i];
            }
            vecI642string(W12_01, temp0);
            message_send[0] = temp0;
#if defined MYDEBUG
            // printf("[LOG]:\tP1发送者 通信前\n");
            // print_plus();
            // printf("C12\n");
            // printVecI64(C12);
            // printf("W12_0\n");
            // printVecI64(W12_0);
            // printf("W12_1\n");
            // printVecI64(W12_1);
            // printf("W12_01\n");
            // printVecI64(W12_01);
            // print_minus();
#endif
        }
        messenger->send_and_recv(message_send, message_recv);
        array<VecI64, 2> A0j, A12j;
        for (size_t i = 0; i < 2; i++)
        {
            A0j[i].resize(vec_size);
            A12j[i].resize(vec_size);
        }
        // P0发送者 通信后
        {
            VecI64 MW01_0(vec_size), MW01_1(vec_size);
            assert(message_recv[0].size() == 3 * vec_size * sizeof(int64_t));
            assert(message_recv[1].size() == 3 * vec_size * sizeof(int64_t));

            temp0 = message_recv[0].substr(0, vec_size * sizeof(int64_t));
            temp1 = message_recv[0].substr(vec_size * sizeof(int64_t), vec_size * sizeof(int64_t));
            string2vecI64(temp0, MW01_0);
            string2vecI64(temp1, MW01_1);
            VecI64 W01_01(vec_size);
            temp0 = message_recv[1].substr(0, vec_size * sizeof(int64_t));
            string2vecI64(temp0, W01_01);
            for (size_t i = 0; i < vec_size; i++)
            {
                A0j[1][i] = C02[i];
                A0j[0][i] = Bj[0][i] == 1 ? MW01_1[i] ^ W01_01[i] : MW01_0[i] ^ W01_01[i];
            }
#if defined MYDEBUG
            // printf("[LOG]:\tP0发送者 通信后\n");
            // print_plus();
            // printf("MW01_0\n");
            // printVecI64(MW01_0);
            // printf("MW01_1\n");
            // printVecI64(MW01_1);
            // printf("W01_01\n");
            // printVecI64(W01_01);
            // printf("A0j[0]\n");
            // printVecI64(A0j[0]);
            // printf("A0j[1]\n");
            // printVecI64(A0j[1]);
            // print_minus();
#endif
        }
        // P1发送者 通信后
        {
            VecI64 MW10_0(vec_size), MW10_1(vec_size);
            assert(message_recv[0].size() == 3 * vec_size * sizeof(int64_t));
            assert(message_recv[1].size() == 3 * vec_size * sizeof(int64_t));
            temp0 = message_recv[1].substr(vec_size * sizeof(int64_t), vec_size * sizeof(int64_t));
            temp1 = message_recv[1].substr(2 * vec_size * sizeof(int64_t), vec_size * sizeof(int64_t));
            string2vecI64(temp0, MW10_0);
            string2vecI64(temp1, MW10_1);
            VecI64 W10_01(vec_size);
            temp0 = message_recv[0].substr(2 * vec_size * sizeof(int64_t), vec_size * sizeof(int64_t));
            string2vecI64(temp0, W10_01);
            for (size_t i = 0; i < vec_size; i++)
            {
                A12j[0][i] = C12[i];
                A12j[1][i] = Bj[1][i] == 1 ? MW10_1[i] ^ W10_01[i] : MW10_0[i] ^ W10_01[i];
            }
#if defined MYDEBUG
            // printf("[LOG]:\tP1发送者 通信后\n");
            // print_plus();
            // printf("MW10_0\n");
            // printVecI64(MW10_0);
            // printf("MW10_1\n");
            // printVecI64(MW10_1);
            // printf("W10_01\n");
            // printVecI64(W10_01);
            // printf("A12j[0]\n");
            // printVecI64(A12j[0]);
            // printf("A12j[1]\n");
            // printVecI64(A12j[1]);
            // print_minus();
#endif
        }
        for (size_t i = 0; i < 2; i++)
        {
            Cj[i] = A0j[i] + A12j[i];
        }
#if defined MYDEBUG
        vecbool B_debug;
        reveal_bools(Bj, B_debug);
        VecI64 A0_debug;
        reveal(A0j, A0_debug);
        if (party_id == 0)
            for (size_t i = 0; i < vec_size; i++)
            {
                assert(A0_debug[i] == Aj[0][i] * B_debug[i]);
            }
        VecI64 A12_debug;
        reveal(A12j, A12_debug);
        if (party_id == 1)
            for (size_t i = 0; i < vec_size; i++)
            {
                assert(A12_debug[i] == (Aj[0][i] + Aj[1][i]) * B_debug[i]);
            }

#endif
    }
}
