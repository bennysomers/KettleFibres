import numpy as np
import IAPWS_IF97 as iapws


def refractive_index_thormahlen(rho, T, lamda):
    # Calculates the refractive index of water using the formula provided in
    # Thormahlen et al (1985)

    # Args: rho - density in kg/m^3
    #       T - temperature in Kelvin
    #       lamda - wavelength of light in nm
    # Returns: n - refractive index of water

    # The coefficients:
    a = np.array([3.036167E-03,
                  5.242100E-02,
                  2.117579E-01,
                  -5.195756E-02,
                  5.922248E-02,
                  -1.918429E-02,
                  2.582351E-03,
                  -2.352054E-04,
                  3.964628E-05,
                  3.336153E-02,
                  -4.008264E-02,
                  8.339681E-03,
                  -1.054741E-02,
                  9.491575E-03])

    # Normalisation constants
    rho_0 = 1000  # kg/m^3
    lamda_NA = 589.0  # nm
    T0 = 273.15  # K

    # Normalise the input parameters
    rhostar = rho/rho_0
    lamdastar = lamda/lamda_NA
    T = T/T0

    Rm = a[0]/(lamdastar**2-a[1]) + a[2] + \
        (a[3] + a[4]*lamdastar + a[5]*lamdastar**2 + a[6]*lamdastar**3 + a[7]*lamdastar**4)*lamdastar**2 + \
        a[8]/rhostar + (a[9] + a[10]*lamdastar + a[11]*lamdastar**2) * \
        lamdastar**2*T + (a[12] + a[13]*lamdastar)*lamdastar*T**2

    n = np.sqrt((1+2*Rm*rho) / (1-Rm*rho))
    return n


def refractive_index_R9_97(P, T, wavelength):
    # Calculates the refractive index of water at the given pressure, temperature and wavelength
    # using the IAPWS R9-97 model

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    #       wavelength - wavelength of light in nm
    # Returns: n - refractive index of water at the given pressure, temperature and wavelength

    a = np.array([0.244257733,
                  9.74634476e-3,
                  -3.73234996e-3,
                  2.68678472e-4,
                  1.58920570e-3,
                  2.45934259e-3,
                  0.900704920,
                  -1.66626219e-2])

    lamda_uv = 0.2292020
    lamda_ir = 5.432937

    rhostar = 1000  # kg/m^3
    Tstar = 273.15  # K
    lamdastar = 0.5893  # um

    rho = iapws.density(P, T)

    Tbar = T / Tstar
    rhobar = rho / rhostar
    lamdabar = wavelength / 1000 / lamdastar  # convert to um first

    rhs = a[0] + a[1]*rhobar + a[2]*Tbar + a[3]*lamdabar**2*Tbar + a[4]/lamdabar**2 + \
        a[5]/(lamdabar**2 - lamda_uv**2) + a[6] / \
        (lamdabar**2 - lamda_ir**2) + a[7]*rhobar**2

    n = np.sqrt((1+2*rhs*rhobar) / (1-rhs * rhobar))
    return n


def main():
    P = 30  # MPa
    T = 700-273.15  # degrees Celsius
    wavelength = 589.32  # nm

    density = iapws.determine_density(P, T)

    Rm = density_refraction(density, T, wavelength)
    refractive_index = calculate_refractive_index(density, Rm)

    print(f"density: {density}")
    print(f"Rm: {Rm}")
    print(f"refractive index: {refractive_index}")


if __name__ == "__main__":
    main()
