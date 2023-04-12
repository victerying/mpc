#include "PRgenerator.hpp"

static inline uint64_t prvhash_core64(uint64_t *const Seed0, uint64_t *const lcg0, uint64_t *const Hash0)
{
    uint64_t Seed = *Seed0;
    uint64_t lcg = *lcg0;
    uint64_t Hash = *Hash0;

    Seed *= lcg * 2 + 1;
    const uint64_t rs = Seed >> 32 | Seed << 32;
    Hash += rs + 0xAAAAAAAAAAAAAAAA;
    lcg += Seed + 0x5555555555555555;
    Seed ^= Hash;
    const uint64_t out = lcg ^ rs;
    *Seed0 = Seed;
    *lcg0 = lcg;
    *Hash0 = Hash;
    return (out);
}

PRgenerator::PRgenerator(uint64_t _init_seed, uint64_t _init_lcg) : seed(_init_seed), lcg(_init_lcg), hash(0), int_buff(BUFF_NUM), bool_buff(BUFF_NUM)
{
    for (size_t i = 0; i < 5; i++)
    {
        prvhash_core64(&this->seed, &this->lcg, &this->hash);
    }
    for (size_t i = 0; i < BUFF_NUM; i++)
    {
        int_buff[i] = prvhash_core64(&this->seed, &this->lcg, &this->hash);
        uint64_t temp = prvhash_core64(&this->seed, &this->lcg, &this->hash);
        bool_buff[i] = temp % 2;
    }
}
int64_t PRgenerator::pop_int64()
{
    return (int64_t)prvhash_core64(&this->seed, &this->lcg, &this->hash);
}
int8_t PRgenerator::pop_bool()
{
    uint64_t temp = prvhash_core64(&this->seed, &this->lcg, &this->hash);
    return (int8_t)(temp % 2);
}