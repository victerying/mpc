#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "../utils/SocketMessenger.hpp"

int main(int argc, char const *argv[])
{
    if (argc != 2)
    {
        printf("usage:\ntest2 <partyid>\n");
        return 1;
    }
    size_t party_id = strtol(argv[1], NULL, 10);
    if (party_id > 2)
    {
        printf("party_id must be in {1,2,3}\n");
        return 1;
    }

    SocketMessenger socketMessenger(party_id);
    VirtuleMessenger *virtuleMessengerPtr = &socketMessenger;
    vector<string> message_send(3);
    vector<string> message_recv(3);
    assert(message_send.size() == 3);
    for (size_t i = 0; i < 3; i++)
    {
        if (i == party_id)
            continue;
        message_send[i].clear();
        message_send[i] += std::to_string(party_id);
        message_send[i] += string(" send ");
        message_send[i] += std::to_string(i);
    }
    virtuleMessengerPtr->send_and_recv(message_send, message_recv);

    for (size_t i = 0; i < 3; i++)
    {
        if (i == party_id)
            continue;
        printf("recv from %lu :\n", i);
        std::cout << message_recv[i] << "\n";
    }
    return 0;
}