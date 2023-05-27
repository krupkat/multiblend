#include "mb/logging.h"

#include <memory>
#include <string>

#include <spdlog/spdlog.h>

namespace multiblend::utils {

namespace {

std::shared_ptr<spdlog::logger> mb_logger_ = nullptr;

}  // namespace

void SetLogger(std::shared_ptr<spdlog::logger> logger) { mb_logger_ = logger; }

void Info(const std::string& msg) {
  if (mb_logger_) {
    mb_logger_->info(msg);
  }
}

void Warn(const std::string& msg) {
  if (mb_logger_) {
    mb_logger_->warn(msg);
  }
}

}  // namespace multiblend::utils
