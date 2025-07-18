// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <boost/asio.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/thread/thread.hpp>
#include <cstdint>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <vector>

#include "BalancerConnection.h"

namespace utl {
class Logger;
}

namespace dst {

namespace ip = asio::ip;

const int kWorkersDiscoveryPeriod = 15;  // time in seconds between retrying to
                                         // find new workers on the network
class Distributed;
class LoadBalancer
{
 public:
  // constructor for accepting connection from client
  LoadBalancer(Distributed* dist,
               asio::io_context& service,
               utl::Logger* logger,
               const char* ip,
               const char* workers_domain,
               unsigned short port = 1234);
  ~LoadBalancer();
  bool addWorker(const std::string& ip, unsigned short port);
  void updateWorker(const ip::address& ip, unsigned short port);
  void getNextWorker(ip::address& ip, unsigned short& port);
  void removeWorker(const ip::address& ip,
                    unsigned short port,
                    bool lock = true);
  void punishWorker(const ip::address& ip, unsigned short port);

 private:
  struct Worker
  {
    ip::address ip;
    unsigned short port;
    unsigned short priority;
    Worker(ip::address ip, unsigned short port, unsigned short priority)
        : ip(ip), port(port), priority(priority)
    {
    }
    bool operator==(const Worker& rhs) const
    {
      return (ip == rhs.ip && port == rhs.port && priority == rhs.priority);
    }
  };
  struct CompareWorker
  {
    bool operator()(Worker const& w1, Worker const& w2)
    {
      return w1.priority > w2.priority;
    }
  };

  Distributed* dist_;
  tcp::acceptor acceptor_;
  asio::io_context* service_;
  utl::Logger* logger_;
  std::priority_queue<Worker, std::vector<Worker>, CompareWorker> workers_;
  std::mutex workers_mutex_;
  std::unique_ptr<asio::thread_pool> pool_;
  std::mutex pool_mutex_;
  uint32_t jobs_;
  std::atomic<bool> alive_ = true;
  boost::thread workers_lookup_thread_;
  std::vector<std::string> broadcastData_;

  void start_accept();
  void handle_accept(const BalancerConnection::Pointer& connection,
                     const boost::system::error_code& err);
  void lookUpWorkers(const char* domain, unsigned short port);
  friend class dst::BalancerConnection;
};
}  // namespace dst
