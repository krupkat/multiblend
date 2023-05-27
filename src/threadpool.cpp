#include "mb/threadpool.h"

#include <thread>

namespace multiblend::mt {

Threadpool* Threadpool::instance_;

/**********************************************************************
 * Constructor (private)
 **********************************************************************/
Threadpool::Threadpool(int threads) {
  n_threads_ = threads > 0 ? (std::min)((unsigned int)threads,
                                        std::thread::hardware_concurrency())
                           : std::thread::hardware_concurrency();
  threads_ = new tp_struct[n_threads_];

  for (int i = 0; i < n_threads_; ++i) {
#ifdef _WIN32
    threads_[i].handle = CreateThread(
        nullptr, 1, (LPTHREAD_START_ROUTINE)Thread, &threads_[i], 0, nullptr);
#else
    pthread_create(&threads_[i].handle, nullptr, TP_Thread, &threads_[i]);
#endif
    threads_[i].main_mutex = &main_mutex_;
    threads_[i].return_mutex = &return_mutex_;
    threads_[i].main_cond = &main_cond_;
    threads_[i].return_cond = &return_cond_;
    threads_[i].free = true;
    threads_[i].stop = false;
    threads_[i].queue = &queue_;
    threads_[i].i = i;
  }
}

/**********************************************************************
 * Destructor
 **********************************************************************/
Threadpool::~Threadpool() {
  int i;

  {
    std::lock_guard<std::mutex> mlock(main_mutex_);
    for (i = 0; i < n_threads_; ++i) {
      threads_[i].stop = true;
    }
  }

  main_cond_.notify_all();
  for (i = 0; i < n_threads_; ++i) {
#ifdef _WIN32
    WaitForSingleObject(threads_[i].handle, INFINITE);
#else
    pthread_join(threads_[i].handle, nullptr);
#endif
  }
}

/**********************************************************************
 * Threads
 **********************************************************************/
#ifdef _WIN32
DWORD WINAPI Threadpool::Thread(void* param) {
#else
void* TP_Thread(void* param) {
#endif

  auto* thread_state = (Threadpool::tp_struct*)param;

  while (true) {
    {
      std::unique_lock<std::mutex> mlock(*thread_state->main_mutex);
      thread_state->main_cond->wait(mlock, [=] {
        return !thread_state->queue->empty() || thread_state->stop;
      });
      if (!thread_state->queue->empty()) {
        thread_state->function = thread_state->queue->front();
        thread_state->queue->pop_front();
        thread_state->free = false;
      }
    }
    if (thread_state->stop) {
      break;
    }

    thread_state->function();

    {
      std::lock_guard<std::mutex> mlock(
          *thread_state->return_mutex);  // necessary
      thread_state->free = true;
    }
    thread_state->return_cond->notify_all();
  }

#ifdef _WIN32
  return 0;
#else
  return nullptr;
#endif
}

/**********************************************************************
 * Wait
 **********************************************************************/
void Threadpool::Wait() {
  if (queue_.empty()) {
    int i;
    for (i = 0; i < n_threads_; ++i) {
      if (!threads_[i].free) {
        break;
      }
    }
    if (i == n_threads_) {
      return;
    }
  }

  {
    std::unique_lock<std::mutex> rlock(return_mutex_);
    return_cond_.wait(rlock, [this] {
      if (!queue_.empty()) {
        return false;
      }
      for (int i = 0; i < n_threads_; ++i) {
        if (!threads_[i].free) {
          return false;
        }
      }
      return true;
    });
  }
}

/**********************************************************************
 * Queue
 **********************************************************************/
void Threadpool::Queue(std::function<void()> function) {
  std::lock_guard<std::mutex> mlock(main_mutex_);  // not sure what this is for
  queue_.push_back(std::move(function));
  main_cond_.notify_one();  // changed from notify_all()
}

}  // namespace multiblend::mt
