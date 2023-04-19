#pragma once
#include "../utils/VirtuleMessenger.hpp"
#include "../utils/PRgenerator.hpp"
#include "../eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <array>
#include <string>
#define KEY01 0xAB01
#define KEY12 0xBC12
#define KEY20 0xCA20
#define KEY0 0xF0
#define KEY1 0xF0 
#define KEY2 0xF0
static const uint16_t P = 67;
static const uint64_t L_1 = -1;
using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
// 第一个元素表示最低位的秘密分享,每一个uint16_t表示在Z67上的比特分享，属于[0,66]
using zpshare = std::array<uint16_t, 64>;
using std::string;
using std::vector;

class SecurennEngine
{
private:
    /* data */
public:
    VirtuleMessenger *messenger;
    size_t party_id;
    vector<PRgenerator> prgs;
    vector<string> message_send;
    vector<string> message_recv;
    SecurennEngine(VirtuleMessenger *_messenger, size_t _party_id);
    SecurennEngine(const SecurennEngine &other) = delete;
    void matrix_mul(const MatrixXl &Xj, const MatrixXl &Yj, MatrixXl &Zj);
    void vector_mul(const vector<int64_t> &Xj, const vector<int64_t> &Yj, vector<int64_t> &Zj);
    void reveal(const vector<int64_t> &Xj, vector<int64_t> &X_plain);
    void reveal_matrix(const MatrixXl &Xj, MatrixXl &X_plain);
    void reveal_vec_bools(const vector<int8_t> &Xj, vector<int8_t> &X_plain);


    void reveal_L_1(const vector<int64_t> &Xj, vector<int64_t> &X_plain);
    void reveal_zp(const vector<zpshare> &Xbits_j, vector<zpshare> &Xbits_plain);

    void share(const vector<int64_t> &X_plain, vector<int64_t> &Xj);
    void privateCompare(const vector<zpshare> &Xbits_j, const vector<int8_t> &beta1_plain, const vector<int64_t> &r_palin, vector<int8_t> &beta2_plain);
    void shareConvert(const vector<int64_t> &Aj, vector<int64_t> &Yj);
    void msb_L_1(const vector<int64_t> &Aj, vector<int64_t> &Bj);
    void msb(const vector<int64_t> &Aj, vector<int64_t> &Bj);


    ~SecurennEngine();
};
