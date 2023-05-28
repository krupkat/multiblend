#pragma once

#include <BS_thread_pool.hpp>

namespace multiblend::mt {

class Threadpool {
 public:
  void Queue(std::function<void()> function) {
    instance_.push_task(std::move(function));
  }

  [[nodiscard]] int GetNThreads() const {
    return instance_.get_thread_count();
  };

  void Wait() { instance_.wait_for_tasks(); }

 private:
  explicit Threadpool(int threads) : instance_{threads} {}
  BS::thread_pool instance_;

  friend Threadpool* GetInstance(int threads);
};

Threadpool* GetInstance(int threads = 0);

}  // namespace multiblend::mt
