#pragma once

#include <future>

#include <BS_thread_pool.hpp>

namespace multiblend::mt {

using MultiFuture = BS::multi_future<void>;

class Threadpool {
 public:
  template <typename... Args>
  [[nodiscard]] std::future<void> Queue(Args&&... args) {
    return instance_.submit(std::forward<Args>(args)...);
  }

  [[nodiscard]] int GetNThreads() const {
    return instance_.get_thread_count();
  };

 private:
  explicit Threadpool(int threads)
      : instance_{static_cast<BS::concurrency_t>(threads)} {}
  BS::thread_pool instance_;

  friend Threadpool* GetInstance(int threads);
};

Threadpool* GetInstance(int threads = 0);

}  // namespace multiblend::mt
