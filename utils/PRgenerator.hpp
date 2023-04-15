#pragma once
#include "stdint.h"
#include "stdlib.h"
#include <vector>
using std::vector;
static const size_t BUFF_NUM=160000;
class PRgenerator
{
private:
    /* data */
public:
    uint64_t seed;
    uint64_t lcg;
    uint64_t hash;
    vector<int64_t> int_buff;
    vector<int8_t> bool_buff;
    PRgenerator() = delete;

    int64_t pop_int64();
    int8_t pop_bool(); 
    
    PRgenerator(uint64_t _init_seed, uint64_t _init_lcg);
    ~PRgenerator() = default;
};
