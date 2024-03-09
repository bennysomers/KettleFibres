import numpy as np

def calculate_coefficients(T):
    # Calculates the temperature dependent coeeficients of the Sellmeier equation
    #
    # Args: T - temperature in degrees Celsius
    # Returns: coefficients - list of coefficients for the Sellmeier equation

    # The following parameters are for SiO2 glass taken from Ghosh et al (1994) Temperature Dependent Sellmeier Coefficients and Chromatic Dispersions for some Optical Fiber Glasses
    m = np.array([6.90754e-6, 2.35835e-5, 5.84758e-7, 5.48368e-7, 0])
    c1 = np.array([1.31552, 7.88404e-1, 1.10199e-2, 0.91316, 1080])

    coefficients = m*T+c1
    return coefficients

def calculate_refractive_index(wavelength, coefficients):
    # Calculates the refractive index of a material using the Sellmeier equation
    #
    # Args: wavelength - wavelength of light in nm
    #       coefficients - list of coefficients for the Sellmeier equation
    # Returns: n - refractive index of the material at the given wavelength

    A = coefficients[0]
    B = coefficients[1]
    C = coefficients[2]
    D = coefficients[3]
    E = coefficients[4]

    wavelength = wavelength/1000 # convert to micrometers
    n = np.sqrt(A +B/(1- C/wavelength**2)+D/(1- E/wavelength**2))

    return n



def main():
    coefficients = calculate_coefficients(471)
    wavelength = 1550
    refractive_index = calculate_refractive_index(wavelength, coefficients)

    print(coefficients)
    print(refractive_index)

if __name__ == "__main__":
    main()
