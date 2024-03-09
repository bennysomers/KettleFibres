import IAPWS_IF97 as iapws
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# print(iapws.saturation_pressure(425.15))
# print(iapws.saturation_pressure(433.15))

# print(iapws.region3_density(75, 650))


def show_plots(nu, rho):
    plt.subplot(2, 1, 1)
    plt.imshow(nu,
               extent=[0, 800, 0, 100],
               aspect='auto',
               cmap='jet')
    plt.title('Specific volume')
    plt.colorbar()

    plt.subplot(2, 1, 2)
    plt.imshow(rho,
               extent=[0, 800, 0, 100],
               aspect='auto',
               cmap='jet')
    plt.title('Density')
    plt.clim(0, 1000)
    plt.colorbar()
    plt.show()


nBins = 1000

maxP = 100
maxT = 800

P = [float(x) * maxP / nBins for x in range(1, nBins)]
T = [float(x) * maxT / nBins + 273.15 for x in range(nBins)]

nu = np.zeros((nBins-1, nBins))
rho = np.zeros((nBins-1, nBins))
for i in tqdm(range(nBins-1)):  # Pressure
    for j in range(nBins):  # Temperature
        nu[nBins - i-2, j] = iapws.specific_volume(P[i], T[j])
        rho[nBins - i-2, j] = iapws.density(P[i], T[j])

show_plots(nu, rho)
