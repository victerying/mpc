#pragma once
#include "../utils/VirtuleMessenger.hpp"
#include "../utils/PRgenerator.hpp"
#include "../eigen-3.4.0/Eigen/Dense"
#include <vector>
#define KEY01 0xAB01
#define KEY12 0xBC12
#define KEY20 0xCA20
#define KEY0 0xF0
#define KEY1 0xF0
#define KEY2 0xF0
using MatrixXl = Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic>;
using std::vector;
using std::string;
#include <string>
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
    void matrix_mul(const MatrixXl &Xj, const MatrixXl &Yj, MatrixXl &Zj);
    void vector_mul(const vector<int64_t> &Xj, const vector<int64_t> &Yj, vector<int64_t> &Zj);
    ~SecurennEngine();
};
