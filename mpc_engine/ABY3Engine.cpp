#include "ABY3Engine.hpp"

ABY3Engine::ABY3Engine(VirtuleMessenger *_messenger, size_t _party_id):messenger(_messenger), party_id(_party_id)
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
