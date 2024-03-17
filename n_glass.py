import numpy as np
import matplotlib.pyplot as plt


def calculate_coefficients(T):
    # Calculates the temperature dependent coeeficients of the Sellmeier equation
    #
    # Args: T - temperature in degrees Celsius
    # Returns: coefficients - list of coefficients for the Sellmeier equation

    # The following parameters are for SiO2 glass taken from Ghosh et al (1994)
    # Temperature Dependent Sellmeier Coefficients and Chromatic Dispersions
    # for some Optical Fiber Glasses
    m = np.array([6.90754e-6, 2.35835e-5, 5.84758e-7, 5.48368e-7, 0])
    c1 = np.array([1.31552, 7.88404e-1, 1.10199e-2, 0.91316, 1080])

    coefficients = m*T+c1
    return coefficients


def refractive_index(wavelength, T):
    # Calculates the refractive index of a material using the Sellmeier equation
    #
    # Args: wavelength - wavelength of light in nm
    #       temperature - temperature in Kelvin
    # Returns: n - refractive index of the material at the given wavelength

    T = T - 273.15  # convert to degrees Celsius

    coefficients = calculate_coefficients(T)

    A = coefficients[0]
    B = coefficients[1]
    C = coefficients[2]
    D = coefficients[3]
    E = coefficients[4]

    wavelength = wavelength/1000  # convert to micrometers
    n = np.sqrt(A + B/(1 - C/wavelength**2)+D/(1 - E/wavelength**2))

    return n


def main():
    T = range(0, 1400)  # degrees Celsius
    wavelength = 1550
    refractive_index = np.zeros(len(T))
    for i in range(len(T)):
        refractive_index[i] = refractive_index(wavelength, T[i])

    plt.plot(T, refractive_index)
    plt.title('Refractive Index of SiO2 Glass at 1550nm')
    plt.xlabel('Temperature (C)')
    plt.ylabel('Refractive Index')
    plt.show()


if __name__ == "__main__":
    main()
