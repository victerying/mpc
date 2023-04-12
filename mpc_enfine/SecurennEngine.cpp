#include "SecurennEngine.hpp"
#include "omp.h"

#include "../utils/util.hpp"

SecurennEngine::SecurennEngine(VirtuleMessenger *_messenger, size_t _party_id) : messenger(_messenger), party_id(_party_id), prgs()
{
    this->message_send.resize(3);
    this->message_recv.resize(3);

    switch (this->party_id)
    {
    case 0:
        prgs.push_back(PRgenerator(0xF0, 0)); // PRG for self
        prgs.push_back(PRgenerator(KEY01, 0));
        prgs.push_back(PRgenerator(KEY20, 0));
        break;
    case 1:
        prgs.push_back(PRgenerator(KEY01, 0));
        prgs.push_back(PRgenerator(0xF1, 0)); // PRG for self
        prgs.push_back(PRgenerator(KEY12, 0));
        break;
    case 2:
        prgs.push_back(PRgenerator(KEY20, 0));
        prgs.push_back(PRgenerator(KEY12, 0));
        prgs.push_back(PRgenerator(0xF2, 0)); // PRG for self
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
    Zj.resize(rows, cols);

    if (party_id == 2)
    {
        MatrixXl A[2];
        MatrixXl B[2];
        MatrixXl C0{};
        MatrixXl C1{};
#pragma omp parallel num_threads(2)
        {
            int thread_id = omp_get_thread_num();
            A[thread_id].resize(rows, lens);
            B[thread_id].resize(lens, cols);
            for (size_t i = 0; i < rows; i++)
            {
                for (size_t k = 0; k < lens; k++)
                {
                    A[thread_id](i, k) = this->prgs[thread_id].pop_int64();
                }
            }
            for (size_t i = 0; i < lens; i++)
            {
                for (size_t k = 0; k < cols; k++)
                {
                    B[thread_id](i, k) = this->prgs[thread_id].pop_int64();
                }
            }
        } // end parallel

        C0.resize(rows, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t k = 0; k < cols; k++)
            {
                C0(i, k) = this->prgs[2].pop_int64();
            }
        }
        MatrixXl C = (A[0] + A[1]) * (B[0] + B[1]);
        C1 = C - C0;
        matrix2string(C0, message_send[0]);
        matrix2string(C1, message_send[1]);
        messenger->send_and_recv(message_send, message_recv);
    }

    else
    { // for j = 0,1  (party_id)
        MatrixXl Aj(rows, lens);
        MatrixXl Bj(lens, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t k = 0; k < lens; k++)
            {
                Aj(i, k) = this->prgs[2].pop_int64();
            }
        }
        for (size_t i = 0; i < lens; i++)
        {
            for (size_t k = 0; k < cols; k++)
            {
                Bj(i, k) = this->prgs[2].pop_int64();
            }
        }
        MatrixXl Xj_sub_Aj = Xj - Aj;
        MatrixXl Yj_sub_Bj = Yj - Bj;
        string temp1, temp2;
        matrix2string(Xj_sub_Aj, temp1);
        matrix2string(Yj_sub_Bj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear(); // 减少带宽
        messenger->send_and_recv(message_send, message_recv);
        temp1 = message_recv[1 - party_id].substr(0, temp1.size());
        temp2 = message_recv[1 - party_id].substr(temp1.size(), temp2.size());
        MatrixXl X_sub_A(rows, lens);
        MatrixXl Y_sub_B(lens, cols);
        string2matrix(temp1, X_sub_A);
        string2matrix(temp2, Y_sub_B);
        X_sub_A = X_sub_A + Xj_sub_Aj;
        Y_sub_B = Y_sub_B + Yj_sub_Bj;
        MatrixXl Cj(rows, cols);
        string2matrix(message_recv[2], Cj);
        Zj = Cj + (X_sub_A * Bj) + (Y_sub_B * Aj) + (Y_sub_B * X_sub_A * (int64_t)party_id);
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
        vector<int64_t> C0{};
        vector<int64_t> C1{};
#pragma omp parallel num_threads(2)
        {
            int thread_id = omp_get_thread_num();
            A[thread_id].resize(vec_size);
            B[thread_id].resize(vec_size);
            for (size_t i = 0; i < vec_size; i++)
            {
                A[thread_id][i] = this->prgs[thread_id].pop_int64();
                B[thread_id][i] = this->prgs[thread_id].pop_int64();
            }

        } // end parallel

        C0.resize(vec_size);
        C1.resize(vec_size);
#pragma omp parallel for
        for (size_t i = 0; i < vec_size; i++)
        {
            C0[i] = this->prgs[2].pop_int64();
            C1[i] = (A[0][i] + A[1][i]) * (B[0][i] + B[1][i]) - C0[i];
        }
        vector2string(C0, message_send[0]);
        vector2string(C1, message_send[1]);
        messenger->send_and_recv(message_send, message_recv);
    }

    else
    { // for j = 0,1  (party_id)
        vector<int64_t> Aj(vec_size);
        vector<int64_t> Bj(vec_size);
        vector<int64_t> Xj_sub_Aj(vec_size);
        vector<int64_t> Yj_sub_Bj(vec_size);
#pragma omp parallel for
        for (size_t i = 0; i < vec_size; i++)
        {
            Aj[i] = this->prgs[2].pop_int64();
            Bj[i] = this->prgs[2].pop_int64();
            Xj_sub_Aj[i] = Xj[i] - Aj[i];
            Yj_sub_Bj[i] = Yj[i] - Bj[i];
        }
        string temp1, temp2;
        vector2string(Xj_sub_Aj, temp1);
        vector2string(Yj_sub_Bj, temp2);
        message_send[1 - party_id] = temp1 + temp2;
        message_send[2].clear(); // 减少带宽
        messenger->send_and_recv(message_send, message_recv);
        temp1 = message_recv[1 - party_id].substr(0, temp1.size());
        temp2 = message_recv[1 - party_id].substr(temp1.size(), temp2.size());
        vector<int64_t> X_sub_A(vec_size);
        vector<int64_t> Y_sub_B(vec_size);
        string2vector(temp1, X_sub_A);
        string2vector(temp2, Y_sub_B);
        vector<int64_t> Cj(vec_size);
        string2vector(message_recv[2], Cj);
#pragma omp parallel for
        for (size_t i = 0; i < vec_size; i++)
        {
            X_sub_A[i] = X_sub_A[i] + Xj_sub_Aj[i];
            Y_sub_B[i] = Y_sub_B[i] + Yj_sub_Bj[i];
            Zj[i] = Cj[i] + (X_sub_A[i] * Bj[i]) + (Y_sub_B[i] * Aj[i]) + (Y_sub_B[i] * X_sub_A[i] * (int64_t)party_id);
        }
    }
}