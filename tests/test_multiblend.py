import unittest
import os
import subprocess

from PIL import Image
from PIL import ImageChops

MULTIBLEND_PATH = os.path.join("..", "Multiblend.exe")
RESULT_PATH = "result.jpg"


class TestStringMethods(unittest.TestCase):
    def test_multiblend_exists(self):
        self.assertTrue(os.path.exists(MULTIBLEND_PATH))

    def test_multiblend(self):
        if os.path.exists(RESULT_PATH):
            os.remove(RESULT_PATH)

        subprocess.run(
            [
                MULTIBLEND_PATH,
                "-o",
                RESULT_PATH,
                "data/img_0.jpg",
                "-366,816",
                "data/img_1.jpg",
                "-670,809",
                "data/img_2.jpg",
                "-163,772",
            ]
        ).check_returncode()

        result = Image.open(RESULT_PATH)
        expected = Image.open("data/expected.jpg")

        self.assertEqual(result.width, expected.width)
        self.assertEqual(result.height, expected.height)

        diff = ImageChops.difference(result, expected)

        self.assertEqual(diff.getbbox(), None)


if __name__ == "__main__":
    unittest.main()
