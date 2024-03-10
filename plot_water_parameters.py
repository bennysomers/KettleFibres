import IAPWS_IF97 as iapws
import n_water
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
    levels = np.linspace(0, 1050, 1000)
    plt.contourf(T, P, rho, levels, cmap='jet', hold=True)
    plt.clim(0, 1000)
    plt.colorbar()
    levels = np.linspace(0, 1050, 50)
    plt.contour(T, P, rho, levels, colors='k',
                linewidths=0.5, linestyles='solid', hold=True)
    plt.title('Density')
    plt.xlabel('Temperature (Celsius)')
    plt.ylabel('Pressure (MPa)')
    plt.show()


def plot_n(n, P, T):
    # plt.imshow(n,
    #            extent=[0, 800, 0, 100],
    #            aspect='auto',
    #            cmap='jet')
    levels = np.linspace(1, max(n.flatten()), 1000)
    plt.contourf(T, P, n, levels, cmap='jet', hold=True)
    plt.clim(1, max(n.flatten()))
    plt.colorbar()
    levels = np.linspace(1, max(n.flatten()), 50)
    plt.contour(T, P, n, levels, colors='k', linewidths=0.5,
                linestyles='solid', hold=True)
    plt.title('Refractive index @ 1550NM')
    plt.xlabel('Temperature (Celsius)')
    plt.ylabel('Pressure (MPa)')
    plt.show()


nBins = 1000

maxP = 100  # MPa
maxT = 800  # degrees Celsius

P = [float(x) * maxP / nBins for x in range(1, nBins+1)]
T = [float(x) * maxT / nBins + 273.15 for x in range(nBins)]
TCelsius = [float(x) * maxT / nBins for x in range(nBins)]

nu = np.zeros((nBins, nBins))
rho = np.zeros((nBins, nBins))
n = np.zeros((nBins, nBins))
for i in tqdm(range(nBins)):  # Pressure
    for j in range(nBins):  # Temperature
        rho[i, j] = iapws.density(P[i], T[j])
        nu[i, j] = 1 / rho[i, j]

        n[i, j] = \
            n_water.refractive_index_thormahlen(rho[i, j], T[j], 1550)

plot_nu_rho(nu, rho, P, TCelsius)

plot_n(n, P, TCelsius)
