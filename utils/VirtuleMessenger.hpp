#pragma once
#include <vector>
#include <string>
using std::vector;
using std::string;
class VirtuleMessenger
{
private:
public:
    VirtuleMessenger() = default;  
    /// @brief one mpc round, every parties send it's message, and receive message from others
    /// @param message_send message need to be sent, len == 3, one of it is empty string (send to self)
    /// @param message_recv ship received message, len == 3 after function return, one of it is empty string (recv from self)
    virtual void send_and_recv(const vector<string>& message_send,  vector<string>& message_recv) = 0;
    virtual ~VirtuleMessenger() = default;
};
