#pragma once

#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace multiblend::mt {

#ifndef _WIN32
void* TP_Thread(void* param);
#endif

class Threadpool {
 public:
  static Threadpool* GetInstance(int threads = 0) {
    if (instance_ == nullptr) {
      instance_ = new Threadpool(threads);
    }
    return instance_;
  }
  void Queue(std::function<void()> function);
  int GetNThreads() const { return n_threads_; };
  void Wait();

  struct tp_struct {
#ifdef _WIN32
    HANDLE handle;
#else
    pthread_t handle;
#endif
    std::function<void()> function;
    bool free;
    bool stop;
    std::mutex* main_mutex;
    std::mutex* return_mutex;
    std::condition_variable* main_cond;
    std::condition_variable* return_cond;
    std::deque<std::function<void()>>* queue;
    int i;
  };

 private:
  static Threadpool* instance_;
  Threadpool(int threads = 0);  // constructor is private
  ~Threadpool();
#ifdef _WIN32
  static DWORD WINAPI Thread(void* param);
#endif
  tp_struct* threads_;
  std::deque<std::function<void()>> queue_;
  int n_threads_;
  std::mutex main_mutex_;
  std::mutex return_mutex_;
  std::condition_variable main_cond_;
  std::condition_variable return_cond_;
};

}  // namespace multiblend::mt
