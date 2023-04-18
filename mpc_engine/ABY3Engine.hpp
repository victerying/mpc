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

using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
using VecI64 = Eigen::Array<int64_t, Eigen::Dynamic, 1>;
using vecbool = vector<int8_t>;
using std::string;
using std::vector;
using std::array;
using std::bitset;
class ABY3Engine
{
private:
public:
    VirtuleMessenger *messenger;
    size_t party_id;
    vector<PRgenerator> prgs;
    vector<string> message_send;
    vector<string> message_recv;
    ABY3Engine(VirtuleMessenger *_messenger, size_t _party_id);
    void share(const VecI64 &X_plain, array<VecI64, 2> &Xj);
    void share_matrix(const MatrixXl &X_plain, array<MatrixXl, 2> &Xj);
    void share_bools(const vecbool &X_plain, array<vecbool, 2> &Xj);

    void reveal(const array<VecI64, 2> &Xj, VecI64 &X_plain);
    void reveal_matrix(const array<MatrixXl, 2> &Xj, MatrixXl &X_plain);
    void reveal_bools(const array<vecbool, 2> &Xj, vecbool &X_plain);

    void reveal_bits(const array<VecI64, 2> &Xj, VecI64 &X_plain);
    
    void matrix_mul(const array<MatrixXl, 2> &Xj, const array<MatrixXl, 2> &Yj, array<MatrixXl, 2> &Zj);
    void vector_mul(const array<VecI64, 2> &Xj, const array<VecI64, 2> &Yj, array<VecI64, 2> &Zj);
    void vecbool_mul(const array<vecbool, 2> &Xj, const array<vecbool, 2> &Yj, array<vecbool, 2> &Zj);
    void add_prepare(const array<VecI64, 2> &Xj, array<VecI64, 2> &Gj, array<VecI64, 2> &Pj);
    void msb(const array<VecI64, 2> &Gj,const array<VecI64, 2> &Pj, array<vecbool, 2> &Zj);




    ABY3Engine(const ABY3Engine &other) = delete;
    ~ABY3Engine();
};
