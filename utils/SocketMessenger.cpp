#include "SocketMessenger.hpp"
void send_up(int send_socket, const void *buf, size_t buf_len)
{
    size_t ptr = 0;
    int ret;
    while (ptr != buf_len)
    {
        ret = send(send_socket, buf + ptr, (buf_len - ptr) < INT32_MAX ? (buf_len - ptr) : INT32_MAX, 0);
        if (ret <= 0)
        {
            perror("send");
            exit(EXIT_FAILURE);
        }
        ptr += ret;
    }
}
void recv_up(int recv_socket, void *buf, size_t buf_len)
{
    size_t ptr = 0;
    int ret;
    while (ptr != buf_len)
    {
        ret = recv(recv_socket, buf + ptr, (buf_len - ptr) < INT32_MAX ? (buf_len - ptr) : INT32_MAX, 0);
        if (ret <= 0)
        {
            perror("recv");
            exit(EXIT_FAILURE);
        }
        ptr += ret;
    }
}
SocketMessenger::SocketMessenger(size_t _party_id) : party_id{_party_id}
{
    for (size_t i = 0; i < 3; i++)
    {
        recv_fd[i] = -1;
        send_fd[i] = -1;
    }

    int opt = 1;
    struct sockaddr_in address;
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(PORT + party_id);
    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        perror("socket failed");
        exit(EXIT_FAILURE);
    }
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)))
    {
        perror("setsockopt");
        exit(EXIT_FAILURE);
    }
    if (bind(server_fd, (struct sockaddr *)&address, sizeof(address)) < 0)
    {
        perror("bind failed");
        exit(EXIT_FAILURE);
    }
    if (listen(server_fd, 2) < 0)
    {
        perror("listen");
        exit(EXIT_FAILURE);
    }

#pragma omp parallel num_threads(2)
#pragma omp sections
    {
#pragma omp section
        {
            for (size_t i = 0; i < 3; i++)
            {
                if (i == party_id)
                    continue;

                struct sockaddr_in serv_addr;
                serv_addr.sin_family = AF_INET;
                serv_addr.sin_port = htons(PORT + i);

                // Convert IPv4 and IPv6 addresses from text to binary
                if (inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr) <= 0)
                {
                    perror("inet_pton");
                    exit(EXIT_FAILURE);
                }
                send_fd[i] = socket(AF_INET, SOCK_STREAM, 0);
                if (send_fd[i] < 0)
                {
                    perror("socket");
                    exit(EXIT_FAILURE);
                }
                int status;
                while (true)
                {
                    status = connect(send_fd[i], (struct sockaddr *)&serv_addr, sizeof(serv_addr));
                    if (status < 0) // 没链接上
                        continue;
                    else
                        break;
                }
                send_up(send_fd[i], (void *)&party_id, sizeof(size_t));
            }
        }
#pragma omp section
        {
            int new_socket;
            size_t remote_id;
            for (size_t i = 0; i < 2; i++)
            {
                new_socket = accept(server_fd, NULL, NULL);
                if (new_socket < 0)
                {
                    perror("accept");
                    exit(EXIT_FAILURE);
                }
                recv_up(new_socket, (void *)&remote_id, sizeof(size_t));
                if (remote_id > 2 || remote_id == party_id)
                {
                    printf("wrong remote_id\n");
                    exit(EXIT_FAILURE);
                }
                recv_fd[remote_id] = new_socket;
            }
        }
    } // 有一个隐藏的barrier
}

SocketMessenger::~SocketMessenger()
{
    for (size_t i = 0; i < 3; i++)
    {
        if (i == party_id)
            continue;
        close(send_fd[i]);
        close(recv_fd[i]);
    }
    shutdown(server_fd, SHUT_RDWR);
}

void SocketMessenger::send_and_recv(const vector<string> &message_send, vector<string> &message_recv)
{
    // send
    for (size_t i = 0; i < 3; i++)
    {
        if (i == party_id)
            continue;
        size_t message_size;
        message_size = message_send[i].size();
        send_up(send_fd[i], (const void *)&message_size, sizeof(size_t));
        send_up(send_fd[i], (const void *)message_send[i].c_str(), message_size);
    }
    // recv
    for (size_t i = 0; i < 3; i++)
    {
        if (i == party_id)
            continue;
        size_t message_size;
        recv_up(recv_fd[i], (void *)&message_size, sizeof(size_t));
        message_recv[i].resize(message_size);
        recv_up(recv_fd[i], (void *)message_recv[i].c_str(), message_size);
    }
}