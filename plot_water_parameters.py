import IAPWS_IF97 as iapws
import n_water
import n_glass
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# print(iapws.saturation_pressure(425.15))
# print(iapws.saturation_pressure(433.15))

# print(iapws.region3_density(75, 650))


def plot_nu_rho(nu, rho, P, T):
    # plt.subplot(2, 1, 1)
    # # # plt.imshow(np.log(nu),
    # #            extent=[0, 800, 0, 100],
    # #            aspect='auto',
    # #            cmap='jet')
    # levels = np.linspace(min(np.log(nu.flatten())), max(np.log(nu.flatten())), 1000)
    # plt.contourf(P, T, np.log(nu), levels, cmap='jet', hold=True)
    # levels = np.linspace(min(np.log(nu.flatten())), max(np.log(nu.flatten())), 50)
    # plt.contour(P, T, np.log(nu), levels, colors='k',
    #             linewidths=0.5, linestyles='solid', hold=True)
    # plt.title('Specific volume (Log scale)')
    # plt.colorbar()

    # plt.subplot(2, 1, 2)
    # plt.imshow(rho,
    #            extent=[0, 800, 0, 100],
    #            aspect='auto',
    #            cmap='jet')

    plt.figure(figsize=(12, 6))
    levels = np.linspace(0, max(rho.flatten()), 1000)
    plt.contourf(T, P, rho, levels, cmap='jet', hold=True)
    plt.clim(0, 1000)
    plt.colorbar(ticks=range(0, 1100, 100), label=r"Density ($kg/m^3$)")
    levels = np.linspace(0, 1050, 22)
    plt.contour(T, P, rho, levels, colors='k',
                linewidths=0.5, linestyles='solid', hold=True)
    plt.title('Density of water')
    plt.xlabel(r'Temperature ($\degree C$)')
    plt.ylabel(r'Pressure ($MPa$)')
    plt.show()


def plot_n(n, P, T):
    # plt.imshow(n,
    #            extent=[0, 800, 0, 100],
    #            aspect='auto',
    #            cmap='jet')
    plt.figure(figsize=(12, 6))
    levels = np.linspace(1, max(n.flatten()), 1000)
    plt.contourf(T, P, n, levels, cmap='jet', hold=True)
    plt.clim(1, max(n.flatten()))
    plt.colorbar(ticks=np.arange(1.0, 1.4, 0.05), label=r"Refractive index")
    levels = np.linspace(1, 1.35, 35)
    plt.contour(T, P, n, levels, colors='k', linewidths=0.5,
                linestyles='solid', hold=True)
    plt.title(r'Refractive index @ $1550nm$')
    plt.xlabel(r'Temperature ($\degree C$)')
    plt.ylabel(r'Pressure ($MPa$)')
    plt.show()


def plot_reflectance(reflectance, P, T):
    reflectance_dB = 10*np.log10(reflectance)
    plt.figure(figsize=(12, 6))
    levels = np.linspace(min(reflectance_dB.flatten()),
                         max(reflectance_dB.flatten()), 1000)
    plt.contourf(T, P, reflectance_dB, cmap='jet', levels=levels, hold=True)
    plt.clim(-28, -14)
    plt.colorbar(ticks=range(-28, -14, 2), label=r"Reflectance ($dB$)")
    levels = np.linspace(-28,
                         -14, 14)
    plt.contour(T, P, reflectance_dB, colors='k', linewidths=0.5,
                linestyles='solid', levels=levels,  hold=True)
    plt.title(r'Reflectance @ $1550nm$')
    plt.xlabel(r'Temperature ($\degree C$)')
    plt.ylabel(r'Pressure ($MPa$)')
    plt.show()


def reflectance(n1, n2):
    return ((n1 - n2) / (n1 + n2))**2


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


nBins = 100

maxP = 100  # MPa
maxT = 800  # degrees Celsius

P = [float(x) * maxP / nBins for x in range(1, nBins+1)]
T = [float(x) * maxT / nBins + 273.15 for x in range(nBins)]
TCelsius = [float(x) * maxT / nBins for x in range(nBins)]

# P = [float(x) * 3/nBins + 19 for x in range(nBins)]
# T = [float(x) * 10/nBins + 370 + 273.15 for x in range(nBins)]
# TCelsius = [float(x) * 3/nBins + 19 + 273.15 for x in range(nBins)]

nu = np.zeros((nBins, nBins))
rho = np.zeros((nBins, nBins))
_n_water = np.zeros((nBins, nBins))
_n_glass = np.zeros(nBins)
_reflectance = np.zeros((nBins, nBins))
for j in tqdm(range(nBins)):  # Temperature
    _n_glass[j] = n_glass.calculate_refractive_index(1550, TCelsius[j])
    for i in range(nBins):  # Pressure
        rho[i, j] = iapws.density(P[i], T[j])
        nu[i, j] = 1 / rho[i, j]

        _n_water[i, j] = \
            n_water.refractive_index_thormahlen(rho[i, j], T[j], 1550)

        _reflectance[i, j] = reflectance(_n_water[i, j], _n_glass[j])

plot_nu_rho(nu, rho, P, TCelsius)

plot_n(_n_water, P, TCelsius)

plot_reflectance(_reflectance, P, TCelsius)

plot_temp_reflectance(_reflectance, TCelsius, P, P_desired=15)

plot_density_reflectance(_reflectance, rho, P, P_desired=15)

plot_specific_volume_reflectance(_reflectance, nu, P, P_desired=15)
