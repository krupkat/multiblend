// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mb/multiblend.h"

#include <fstream>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include "mb/image.h"
#include "mb/threadpool.h"

multiblend::io::Image Load(const char* path) {
  std::ifstream ifs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(ifs);

  multiblend::io::InMemoryImage image;
  archive(image);
  return multiblend::io::Image(image);
}

TEST_CASE("MultiblendLibTest") {
  std::vector<multiblend::io::Image> images;
  images.push_back(Load("data/img_0.png.raw"));
  images.push_back(Load("data/img_1.png.raw"));
  images.push_back(Load("data/img_2.png.raw"));

  multiblend::Options options{
      .output_type = multiblend::io::ImageType::MB_IN_MEMORY, .output_bpp = 8};

  auto threadpool = multiblend::mt::Threadpool{};

  auto result = multiblend::Multiblend(
      images, options, multiblend::mt::ThreadpoolPtr{&threadpool});

  CHECK(result.width == 331);
  CHECK(result.height == 244);
  CHECK(result.output_bpp == 8);
  CHECK(result.no_mask == false);
  CHECK(result.min_xpos == -167);
  CHECK(result.min_ypos == 193);
}
