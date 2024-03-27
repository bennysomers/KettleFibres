import IAPWS_IF97
import matplotlib.pyplot as plt
import numpy as np
import n_glass
import n_water


def plot_temp_reflectance(_reflectance, TCelsius, P, P_desired=15):
    # Plots the relationship between temperature and reflectance at a desired pressure

    # Find the index of the desired pressure
    i = 0
    for j in range(len(P)):
        if P[j] > P_desired:
            i = j
            break

    reflectance_dB = 10*np.log10(_reflectance[i, :])
    plt.figure(figsize=(12, 4))
    plt.plot(reflectance_dB, TCelsius, 'x-')
    plt.title(f'Temperature @ {P_desired}MPa, $1550nm$')
    plt.ylabel(r'Temperature ($\degree C$)')
    plt.xlabel(r'Reflectance ($dB$)')
    plt.xlim(min(reflectance_dB), max(reflectance_dB)+0.5)
    plt.ylim(min(TCelsius), max(TCelsius))
    plt.show()


def plot_density_reflectance(_reflectance, density, P, P_desired=15):
    # Plots the relationship between temperature and reflectance at a desired pressure

    # Find the index of the desired pressure
    i = 0
    for j in range(len(P)):
        if P[j] > P_desired:
            i = j
            break

    reflectance_dB = 10*np.log10(_reflectance[i, :])
    density = density[i, :]
    plt.figure(figsize=(12, 4))
    plt.plot(reflectance_dB, density, 'x-')
    plt.title(f'Density @ {P_desired}MPa, $1550nm$')
    plt.ylabel(r'Density ($kg/m^3$)')
    plt.xlabel(r'Reflectance ($dB$)')
    plt.xlim(min(reflectance_dB), max(reflectance_dB)+0.5)
    plt.ylim(min(density), max(density))
    plt.show()


def plot_specific_volume_reflectance(_reflectance, nu, P, P_desired=15):
    # Plots the relationship between temperature and reflectance at a desired pressure

    # Find the index of the desired pressure
    i = 0
    for j in range(len(P)):
        if P[j] > P_desired:
            i = j
            break

    reflectance_dB = 10*np.log10(_reflectance[i, :])
    nu = nu[i, :]
    plt.figure(figsize=(12, 4))
    plt.plot(reflectance_dB, nu, 'x-')
    plt.title(f'Specific Volume @ {P_desired}MPa, $1550nm$')
    plt.ylabel(r'Density ($m^3/kg$)')
    plt.xlabel(r'Reflectance ($dB$)')
    plt.xlim(min(reflectance_dB), max(reflectance_dB)+0.5)
    plt.ylim(min(nu), max(nu))
    plt.show()


R = [12 * x - 27 for x in range(100)]
T = np.zeros(100)
rho = np.zeros(100)
nu = np.zeros(100)


def reflectance(n1, n2):
    return ((n1 - n2) / (n1 + n2))**2


rho = np.zeros(len(R))
nu = np.zeros(len(R))
T = np.zeros(len(R))

for i in range(len(R)):
    # Assume a constant value for pressure
    P = 15  # MPa
    # Guess a value for temperature
    T[i] = 400 + 273.15  # Kelvin

    def temperature_reflectance_residual(_T):
        # calculate the density for the given temperature and pressure
        _rho = IAPWS_IF97.density(P, _T)
        # calculate the refractive index for the given temperature and density
        _n_water = n_water.refractive_index_thormahlen(P, _T, 1550)
        _n_glass = n_glass.refractive_index(1550, _T-273.15)
        # calculate the reflectance for the given refractive indices
        R_calc = reflectance(_n_glass, _n_water)
        # calculate the residual between the calculated reflectance and the measured reflectance
        return R[i] - R_calc

    # Use newtons method to solve for the temperature
    T[i] = IAPWS_IF97.newtons_method(temperature_reflectance_residual, T[i])
    # calculate the density for the given temperature and pressure
    rho[i] = IAPWS_IF97.density(P, T[i])
    nu[i] = 1/rho[i]

plot_specific_volume_reflectance(R, nu, P)
plot_density_reflectance(R, rho, P)
plot_temp_reflectance(R, T, P)
