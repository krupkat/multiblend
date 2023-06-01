// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <future>

#include <BS_thread_pool.hpp>

namespace multiblend::mt {

template <typename TReturnType>
using MultiFuture = BS::multi_future<TReturnType>;

class ThreadpoolPtr;

class ThreadpoolApiWrapper {
 public:
  template <typename... Args>
  [[nodiscard]] auto Queue(Args&&... args) {
    return instance_->submit(std::forward<Args>(args)...);
  }

  [[nodiscard]] int GetNThreads() const {
    return instance_->get_thread_count();
  };

 private:
  explicit ThreadpoolApiWrapper(BS::thread_pool* instance)
      : instance_{instance} {}
  BS::thread_pool* instance_;

  friend ThreadpoolPtr;
};

class ThreadpoolPtr {
 public:
  explicit ThreadpoolPtr(BS::thread_pool* instance) : instance_{instance} {}
  ThreadpoolApiWrapper* operator->() { return &instance_; }

 private:
  ThreadpoolApiWrapper instance_;
};

using Threadpool = BS::thread_pool;

}  // namespace multiblend::mt
