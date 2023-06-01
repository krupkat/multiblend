// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

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

void SetLogger(std::shared_ptr<spdlog::logger> logger) {
  mb_logger_ = std::move(logger);
}

void SetVerbosity(spdlog::logger* logger, int verbosity) {
  switch (verbosity) {
    case 0:
      logger->set_level(spdlog::level::warn);
      return;
    case 1:
      logger->set_level(spdlog::level::info);
      return;
    case 2:
      logger->set_level(spdlog::level::debug);
      return;
    default:
      break;
  }
  if (verbosity < 0) {
    logger->set_level(spdlog::level::err);
  }
  if (verbosity > 2) {
    logger->set_level(spdlog::level::trace);
  }
}

}  // namespace multiblend::utils
