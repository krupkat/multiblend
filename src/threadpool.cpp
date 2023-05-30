#include "mb/threadpool.h"

#include <thread>

namespace multiblend::mt {

Threadpool* GetInstance(int threads) {
  static Threadpool instance{threads};
  return &instance;
}

}  // namespace multiblend::mt
