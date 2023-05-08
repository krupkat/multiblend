#include "src/threadpool.h"

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
    threads_[i].handle = CreateThread(NULL, 1, (LPTHREAD_START_ROUTINE)Thread,
                                      &threads_[i], 0, NULL);
#else
    pthread_create(&threads_[i].handle, NULL, TP_Thread, &threads_[i]);
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
    pthread_join(threads_[i].handle, NULL);
#endif
  }
}

/**********************************************************************
 * Threads
 **********************************************************************/
#define P ((Threadpool::tp_struct*)param)

#ifdef _WIN32
DWORD WINAPI Threadpool::Thread(void* param) {
#else
void* TP_Thread(void* param) {
#endif
  while (true) {
    {
      std::unique_lock<std::mutex> mlock(*P->main_mutex);
      P->main_cond->wait(mlock, [=] { return P->queue->size() || P->stop; });
      if (P->queue->size()) {
        P->function = P->queue->front();
        P->queue->pop_front();
        P->free = false;
      }
    }
    if (P->stop) {
      break;
    }

    P->function();

    {
      std::lock_guard<std::mutex> mlock(*P->return_mutex);  // necessary
      P->free = true;
    }
    P->return_cond->notify_all();
  }

  return 0;
}

/**********************************************************************
 * Wait
 **********************************************************************/
void Threadpool::Wait() {
  if (queue_.size() == 0u) {
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
      if (queue_.size() != 0u) {
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
