import unittest

from gwskysim.detectors.gwdetector import GWDetector


class TestCases(unittest.TestCase):

    def test_init(self):
        # random information picking
        self.assertEqual(GWDetector('V1').longitude, -10.5)

        # it is possible to create every detector in the available_detectors dictionary
        for det_name in GWDetector.get_detector_names():
            GWDetector(det_name)

        # if name is spelled wrong raise an error
        self.assertRaises(ValueError, GWDetector, 'detector_that_does_not_exist')

    def test_get_detector_names(self):
        self.assertEqual(GWDetector.get_detector_names(), GWDetector.available_detectors.keys())

    def test_new_detector(self):
        # creating instance
        thomas = GWDetector.new_detector('thomas_the_detector', 0, 0, 0, 90)

        # thomas is a new detector
        self.assertTrue('thomas_the_detector' in GWDetector.get_detector_names())

        # random info picking
        self.assertEqual(GWDetector('thomas_the_detector').longitude, 0)

        # it is possible to create every detector instance
        for det_name in GWDetector.get_detector_names():
            GWDetector(det_name)


if __name__ == '__main__':
    unittest.main()
