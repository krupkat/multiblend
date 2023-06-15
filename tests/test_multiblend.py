# SPDX-FileCopyrightText: 2023 Tomas Krupka
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest
import os
import subprocess

from PIL import Image
from PIL import ImageChops

MULTIBLEND_PATH = os.path.join("..", "Multiblend")


class TestMultiblend(unittest.TestCase):
    def assert_equal(self, result_path, expected_path):
        with Image.open(result_path) as result, Image.open(expected_path) as expected:
            self.assertEqual(result.width, expected.width)
            self.assertEqual(result.height, expected.height)

            diff = ImageChops.difference(result, expected)
            self.assertEqual(diff.getbbox(), None)

    def test_multiblend_jpg(self):
        result_path = "result.jpg"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "-o",
                result_path,
                "data/img_0.jpg",
                "-366,816",
                "data/img_1.jpg",
                "-670,809",
                "data/img_2.jpg",
                "-163,772",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected.jpg")

    def test_multiblend_png(self):
        result_path = "result.png"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "-o",
                result_path,
                "data/img_0.png",
                "-91,204",
                "data/img_1.png",
                "-167,202",
                "data/img_2.png",
                "-41,193",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected.png")

    def test_multiblend_png_alpha(self):
        result_path = "result_alpha.png"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "-o",
                result_path,
                "data/img_0_alpha.png",
                "-91,204",
                "data/img_1_alpha.png",
                "-167,202",
                "data/img_2_alpha.png",
                "-41,193",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected_alpha.png")

    def test_multiblend_png_saveseams(self):
        result_path = "result_saveseams.png"
        seams_path = "seams.png"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "--saveseams",
                seams_path,
                "-o",
                result_path,
                "data/img_0.png",
                "-91,204",
                "data/img_1.png",
                "-167,202",
                "data/img_2.png",
                "-41,193",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected.png")
        self.assert_equal(seams_path, "data/expected_seams.png")

    def test_multiblend_png_loadseams(self):
        result_path = "result_loadseams.png"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "--loadseams",
                "data/expected_seams.png",
                "-o",
                result_path,
                "data/img_0.png",
                "-91,204",
                "data/img_1.png",
                "-167,202",
                "data/img_2.png",
                "-41,193",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected.png")

    def test_multiblend_tif_8bit(self):
        result_path = "result.tif"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "-o",
                result_path,
                "data/img_0.tif",
                "76,11",
                "data/img_1.tif",
                "0,9",
                "data/img_2.tif",
                "126,0",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected.tif")

    def test_multiblend_tif_16bit(self):
        result_path = "result_16bit.tif"

        if os.path.exists(result_path):
            os.remove(result_path)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "-o",
                result_path,
                "data/img_0_16bit.tif",
                "76,11",
                "data/img_1_16bit.tif",
                "0,9",
                "data/img_2_16bit.tif",
                "126,0",
            ]
        ).check_returncode()

        self.assert_equal(result_path, "data/expected_16bit.tif")


if __name__ == "__main__":
    unittest.main()
