import unittest
from IAPWS_R6_95 import *


class TestIAPWS(unittest.TestCase):

    def test_helmholtz_free_energy_phi_0(self):
        # Test the helmholtz_free_energy_phi_0 function
        T = 500  # K
        rho = 838.025  # kg/m^3
        expected_phi_0 = 0.204797733e1  # kJ/kg
        phi_0 = helmholtz_free_energy_phi_0(T, rho)
        self.assertAlmostEqual(phi_0, expected_phi_0, places=6)

    def test_helmholtz_free_energy_phi_r(self):
        # Test the helmholtz_free_energy_phi_r function
        T = 500  # K
        rho = 838.025  # kg/m^3
        expected_phi_r = -0.342693206e1  # kJ/kg
        phi_r = helmholtz_free_energy_phi_r(T, rho)
        self.assertAlmostEqual(phi_r, expected_phi_r, places=6)

    def test_helmholtz_free_energy_phi_r_derivative_delta(self):
        # Test the helmholtz_free_energy_phi_r_derivative_delta function
        T = 500  # K
        rho = 838.025  # kg/m^3
        expected_derivative = -0.364366650  # kJ/(kg m^3)
        derivative = helmholtz_free_energy_phi_r_derivative_delta(T, rho)
        self.assertAlmostEqual(derivative, expected_derivative, places=6)

    def test_pressure(self):
        # Test the pressure function

        # Test case 1
        T = 300  # K
        rho = 0.996556e3  # kg/m^3
        expected_P = 0.992418352e-1  # MPa
        P = pressure(T, rho)
        self.assertAlmostEqual(P, expected_P, places=6)

    def test_density(self):
        # Test the density function
        P = 0.992418352e-1  # MPa
        T = 300  # K
        expected_rho = 0.9965560e3  # kg/m^3
        rho = density(P, T)
        self.assertAlmostEqual(rho, expected_rho, places=6)


if __name__ == '__main__':
    unittest.main()
