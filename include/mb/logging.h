// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

namespace multiblend::utils {

void Info(const std::string& msg);
void Warn(const std::string& msg);

template <typename... Args>
void Output(int level, fmt::format_string<Args...> fmt, Args&&... args) {
  if (level == 1) {
    Info(fmt::format(fmt, std::forward<Args>(args)...));
  } else {
    Warn(fmt::format(fmt, std::forward<Args>(args)...));
  }
}

void SetLogger(std::shared_ptr<spdlog::logger> logger);

void SetVerbosity(spdlog::logger* logger, int verbosity);

template <typename... Args>
void Throw(fmt::format_string<Args...> fmt, Args&&... args) {
  throw std::runtime_error(fmt::format(fmt, std::forward<Args>(args)...));
}

}  // namespace multiblend::utils
