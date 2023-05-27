#include "mb/logging.h"

#include <memory>
#include <string>

#include <spdlog/spdlog.h>

namespace multiblend::utils {

namespace {

std::shared_ptr<spdlog::logger> mb_logger_ = nullptr;

}  // namespace

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

void SetLogger(std::shared_ptr<spdlog::logger> logger) { mb_logger_ = logger; }

void SetVerbosity(std::shared_ptr<spdlog::logger>, int verbosity) {
  switch (verbosity) {
    case 0:
      spdlog::set_level(spdlog::level::warn);
      return;
    case 1:
      spdlog::set_level(spdlog::level::info);
      return;
    case 2:
      spdlog::set_level(spdlog::level::debug);
      return;
    default:
      break;
  }
  if (verbosity < 0) {
    spdlog::set_level(spdlog::level::err);
  }
  if (verbosity > 2) {
    spdlog::set_level(spdlog::level::trace);
  }
}

}  // namespace multiblend::utils
