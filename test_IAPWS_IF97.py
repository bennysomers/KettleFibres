import unittest
from IAPWS_IF97 import specific_volume, region3_pressure


class TestDensityWater(unittest.TestCase):
    def test_region1_specific_volume(self):
        # Test case 1: P = 3 MPa, T = 300K
        P1 = 3
        T1 = 300
        expected_result1 = 0.100215168e-2
        result1 = specific_volume(P1, T1)
        self.assertAlmostEqual(result1, expected_result1, places=9)

        # Test case 2: P = 80 MPa, T = 300K
        P2 = 80
        T2 = 300
        expected_result2 = 0.971180894e-3
        result2 = specific_volume(P2, T2)
        self.assertAlmostEqual(result2, expected_result2, places=9)

        # Test case 3: P = 3 MPa, T = 500K
        P3 = 3
        T3 = 500
        expected_result3 = 0.120241800e-2
        result3 = specific_volume(P3, T3)
        self.assertAlmostEqual(result3, expected_result3, places=9)

    def test_region2_specific_volume(self):
        # Test case 1: P = 0.0035 MPa, T = 300K
        P1 = 0.0035
        T1 = 300
        expected_result1 = 0.394913866e2
        result1 = specific_volume(P1, T1)
        self.assertAlmostEqual(result1, expected_result1, places=7)

        # Test case 2: P = 0.0035 MPa, T = 700K
        P2 = 0.0035
        T2 = 700
        expected_result2 = 0.923015898e2
        result2 = specific_volume(P2, T2)
        self.assertAlmostEqual(result2, expected_result2, places=7)

        # Test case 3: P = 30 MPa, T = 700K
        P3 = 30
        T3 = 700
        expected_result3 = 0.542946619e-2
        result3 = specific_volume(P3, T3)
        self.assertAlmostEqual(result3, expected_result3, places=9)

    def test_region3_pressure(self):
        # Test case 1: T = 650K, rho = 500 kg/m^3
        T1 = 650
        rho1 = 500
        expected_result1 = 0.255837018e2
        result1 = region3_pressure(rho1, T1)
        self.assertAlmostEqual(result1, expected_result1, places=7)

        # Test case 2: T = 650K, rho = 200 kg/m^3
        T2 = 650
        rho2 = 200
        expected_result2 = 0.222930643e2
        result2 = region3_pressure(rho2, T2)
        self.assertAlmostEqual(result2, expected_result2, places=7)

        # Test case 3: T = 750K, rho = 500 kg/m^3
        T3 = 750
        rho3 = 500
        expected_result3 = 0.783095639e2
        result3 = region3_pressure(rho3, T3)
        self.assertAlmostEqual(result3, expected_result3, places=7)


if __name__ == '__main__':
    unittest.main()
