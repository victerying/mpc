#pragma once
#include "VirtuleMessenger.hpp"
#include "stdlib.h"
#include <string.h>
#include <sys/socket.h>
#include <unistd.h>
#include <netinet/in.h>
#include <omp.h>
#include <arpa/inet.h>
#include <assert.h>
// 第i方的服务器socket占用的端口 =PORT + i
#define PORT 8000

class SocketMessenger : public VirtuleMessenger
{
private:
public:
    size_t party_id;
    int server_fd;
    int recv_fd[3];
    int send_fd[3];
    SocketMessenger(size_t _party_id);
    SocketMessenger(const SocketMessenger &others) = delete;
    virtual void send_and_recv(const vector<string> &message_send, vector<string> &message_recv) override;
    virtual ~SocketMessenger();
};
