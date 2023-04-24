#pragma once
#include "../utils/VirtuleMessenger.hpp"
#include "../utils/PRgenerator.hpp"
#include "../eigen-3.4.0/Eigen/Dense"
#include "../utils/util.hpp"
#include <vector>
#include <array>
#include <bitset>

#define KEY01 0xAB01
#define KEY12 0xBC12
#define KEY20 0xCA20
#define KEY0 0xF0
#define KEY1 0xF0
#define KEY2 0xF0

// using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
// using VecI64 = Eigen::Array<int64_t, Eigen::Dynamic, 1>;
// using vecbool = vector<int8_t>;
using std::array;
using std::bitset;
using std::string;
using std::vector;
class MyEngine
{
private:
    /* data */
public:
    VirtuleMessenger *messenger;
    size_t party_id;
    vector<PRgenerator> prgs;
    vector<string> message_send;
    vector<string> message_recv;
    MyEngine(VirtuleMessenger *_messenger, size_t _party_id);
    ~MyEngine();
    void share(const VecI64 &X_plain, VecI64 &Xj);
    void share_matrix(const MatrixXl &X_plain, MatrixXl &Xj);
    void reveal(const VecI64 &Xj, VecI64 &X_plain);
    void reveal_matrix(const MatrixXl &Xj, MatrixXl &X_plain);
    void matrix_mul(const MatrixXl &Xj, const MatrixXl &Yj, MatrixXl &Zj);
    void vector_mul(const VecI64 &Xj, const VecI64 &Yj, VecI64 &Zj);
    void privateCompare(const vector<zpshare> &Xbits_j, const vecbool &beta1_plain, const VecI64 &r_palin, vecbool &beta2_plain);
    void msb(const VecI64 &Aj, vecbool &Bj);
    void bit_inject(const VecI64 &Xj, const vecbool &Aj, VecI64 &Zj);
    void reveal_vec_bools(const vecbool &Xj, vecbool &X_plain);
    void reveal_L_1(const VecI64 &Xj, VecI64 &X_plain);
    void reveal_zp(const vector<zpshare> &Xbits_j, vector<zpshare> &Xbits_plain);
};
