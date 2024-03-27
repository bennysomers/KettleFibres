# This file attempts to calculate the transmutation of isotopes within a glass fibre.
# The following equation is the basis of everything:
#       dN = -lambda*N*dt + lambda_p*N_p*dt - sigma*N*phi*dt + sigma_p*N_p*phi*dt
#       where:
#           N is the number of atoms of the element
#           dN is the change in number of atoms of a particular element
#           lambda is the decay constant of the element
#           sigma is the cross section of the element
#
#           lambda_p is the decay constant of the parent element
#           N_p is the number of atoms of the parent element
#           sigma_p is the cross section of the parent element
#
#           phi is the neutron flux
#           dt is the time step
#
# There are 2 source terms and 2 sink terms for each isotope present.

import numpy as np
import matplotlib.pyplot as plt
import requests
import csv
from os.path import exists
from tqdm import tqdm
import multiprocessing as mp


class Isotope:
    def __init__(self, Z, A, N):
        self.Z = Z
        self.A = A
        # fraction of atoms in material >>>>>>>>>> NOT NUMBER OF NEUTRONS <<<<<<<<<<
        self.N = np.empty((total_time_steps+1))
        self.N[:] = np.nan
        self.N[range(time_step+1)] = 0
        self.N[time_step] = N
        self.N[time_step+1] = N
        self.elementSymbol = f"{element_symbols()[Z-1]}"
        self.elementName = f"{element_names()[Z-1]}"
        self.name = f"{A}{self.elementSymbol}"
        self.key = f"{Z}_{A}"
        self.get_parameters()

    def __str__(self):
        str = f"Element: {self.elementName}\nIsotope: {self.name}\nZ: {self.Z}\nA: {self.A}"
        if hasattr(self, "half_life"):
            str += f"\nHalf-life: {self.half_life} "
        if hasattr(self, "a_fraction"):
            str += f"\n\N{GREEK SMALL LETTER ALPHA} fraction: {self.a_fraction}"
        if hasattr(self, "bp_fraction"):
            str += f"\n\N{GREEK SMALL LETTER BETA}+ fraction: {self.bp_fraction}"
        if hasattr(self, "bm_fraction"):
            str += f"\n\N{GREEK SMALL LETTER BETA}- fraction: {self.bm_fraction}"
        if hasattr(self, "n_xs"):
            str += f"\n(n, ) cross section: {self.n_xs} cm^2"
        if hasattr(self, "n_gamma_xs"):
            str += f"\n(n, \N{GREEK SMALL LETTER GAMMA}) cross section: {self.n_gamma_xs} cm^2"
        if hasattr(self, "n_alpha_xs"):
            str += f"\n(n, \N{GREEK SMALL LETTER ALPHA}) cross section: {self.n_alpha_xs} cm^2"
        if hasattr(self, "n_proton_xs"):
            str += f"\n(n, p) cross section: {self.n_proton_xs} cm^2"
        return str

    def get_parameters(self):
        url = "https://nds.iaea.org/relnsd/v1/data?"
        rad_types = ["a", "bp", "bm", "g"]
        for rad_type in rad_types:
            csv_name = f"decay_parameters/{self.name}_{rad_type}.csv"
            if not exists(csv_name):
                print(f"Downloading data for {self.name} {rad_type} from IAEA...")
                fields = f"fields=decay_rads&nuclides={self.name}&rad_types={rad_type}"
                request_str = f"{url}{fields}"
                response = requests.get(request_str)
                with open(csv_name, "w") as file:
                    writer = csv.writer(file)
                    for line in response.iter_lines():
                        writer.writerow(line.decode('utf-8').split(","))
            self._get_params_from_csv(csv_name, rad_type)

        # (n, *) cross section
        csv_name = f"capture_parameters/thermal_resonance_ranges.csv"
        with open(csv_name, "r") as file:
            reader = csv.reader(file)
            for row in reader:
                if row[0] == 'Element':
                    continue
                if int(row[1]) == self.A and int(row[2]) == self.Z:
                    self.n_xs = 0
                    # Maxw.(n, gamma) (treating as (n,gamma))
                    if row[5] != "":
                        self.n_gamma_xs = float(row[5])/1e24
                        self.n_xs += self.n_gamma_xs
                    # (n, gamma)
                    if row[7] != "":
                        self.n_gamma_xs = float(row[7])/1e24
                        self.n_xs += self.n_gamma_xs
                    # (n, alpha)
                    if row[11] != "":
                        self.n_alpha_xs = float(row[11])/1e24
                        self.n_xs += self.n_alpha_xs
                    # (n, proton)
                    if row[13] != "":
                        self.n_proton_xs = float(row[13])/1e24
                        self.n_xs += self.n_proton_xs
                    if self.n_xs == 0:
                        del (self.n_xs)
                    break

    def _get_params_from_csv(self, file, rad_type):
        with open(file, "r") as file:
            reader = csv.reader(file)
            line_count = 0
            for row in reader:
                if row[0] == 0:
                    break
                if line_count == 1:
                    if rad_type == "bm":
                        self.half_life = float(row[23])
                        self.decay_constant = np.log(2)/self.half_life
                        try:
                            self.bm_fraction = float(row[26])/100
                        except:
                            self.bm_fraction = 1
                        break
                    elif rad_type == "g":
                        break
                    elif rad_type == "bp":
                        self.half_life = float(row[25])
                        self.decay_constant = np.log(2)/self.half_life
                        try:
                            self.bp_fraction = float(row[28])/100
                        except:
                            self.bp_fraction = 1
                        break
                    elif rad_type == "a":
                        self.half_life = float(row[18])
                        self.decay_constant = np.log(2)/self.half_life
                        try:
                            self.a_fraction = float(row[21])/100
                        except:
                            self.a_fraction = 1
                        break
                line_count += 1

    def decay(self, dt):
        if hasattr(self, "decay_constant"):
            # print(f"============================\n{self.name} decaying...")
            dN = -self.decay_constant*self.N[time_step]*dt
            if dN < -self.N[time_step]:
                dN = -self.N[time_step]
            # print(f"{self.name} has {self.N[time_step]} atoms")
            # print(f"{self.name} decayed by {dN} atoms")
            self.N[time_step+1] += dN
            # print(f"{self.name} has {self.N[time_step+1]} atoms")
            self._decay_to(-dN)

    def _decay_to(self, dN):
        global topes
        # ( , a) decay
        if hasattr(self, "a_fraction"):
            daughter_Z = self.Z - 2
            daughter_A = self.A - 4
            self._add_to_isotope(dN*self.a_fraction, daughter_Z, daughter_A)

            alpha_Z = 2
            alpha_A = 4
            self._add_to_isotope(dN*self.a_fraction, alpha_Z, alpha_A)

        # (, bp) decay
        if hasattr(self, "bp_fraction"):
            daughter_Z = self.Z - 1
            daughter_A = self.A
            self._add_to_isotope(dN*self.bp_fraction, daughter_Z, daughter_A)

        # (, bm) decay
        if hasattr(self, "bm_fraction"):
            daughter_Z = self.Z + 1
            daughter_A = self.A
            self._add_to_isotope(dN*self.bm_fraction, daughter_Z, daughter_A)

    def neutron_capture(self, dt, phi, path_length):
        # This assumes that neutron flux is unidirectional and passes through the full thickness of the fibre.
        # This is a simplification, but it is a good starting point.

        # (n, gamma) capture
        if hasattr(self, "n_xs"):
            # print(f"============================\n{self.name} capturing neutrons...")
            # print(f"{self.name} capturing neutrons...")
            dN = -self.n_xs*self.N[time_step]*phi*dt  # *path_length
            if dN < -self.N[time_step]:
                dN = -self.N[time_step]
            # print(f"{self.name} has {self.N[time_step]} atoms")
            # print(f"{self.name}.n_xs = {self.n_xs}")
            # print(f"{self.name} captured {-dN} neutrons")
            self.N[time_step+1] += dN
            # print(f"{self.name} has {self.N[time_step+1]} atoms")
            self._n_capture_to(-dN)

    def _n_capture_to(self, dN):
        global topes

        ###############################
        # (n, gamma) capture
        if hasattr(self, "n_gamma_xs"):
            daughter_Z = self.Z
            daughter_A = self.A + 1
            dN_gamma = dN*self.n_gamma_xs/self.n_xs
            self._add_to_isotope(dN_gamma, daughter_Z, daughter_A)

        ###############################
        # (n, alpha) capture
        if hasattr(self, "n_alpha_xs"):
            daughter_Z = self.Z-2
            daughter_A = self.A - 4
            dN_alpha = dN*self.n_alpha_xs/self.n_xs
            self._add_to_isotope(dN_alpha, daughter_Z, daughter_A)
            self._add_to_isotope(dN_alpha, 2, 4)
        ###############################
        # (n, proton) capture
        if hasattr(self, "n_proton_xs"):
            daughter_Z = self.Z + 1
            daughter_A = self.A + 1
            dN_proton = dN*self.n_proton_xs/self.n_xs
            self._add_to_isotope(dN_proton, daughter_Z, daughter_A)

    def _add_to_isotope(self, dN, Z, A):
        global topes
        key = f"{Z}_{A}"
        if key not in topes:
            # print(f"adding {key} to topes")
            topes[key] = Isotope(Z, A, 0)
        #     print(f"{topes[key].name} has been added to topes")
        # print(f"{topes[key].name} has {topes[key].N[time_step]} atoms")
        # print(f"{dN} {self.name} atoms were transmuted to {topes[key].name}")
        # if np.isnan(topes[key].N[time_step+1]):
        #     topes[key].N[time_step+1] = topes[key].N[time_step]
        topes[key].N[time_step+1] += dN
        # print(f"{topes[key].name} has {topes[key].N[time_step+1]} atoms")


def element_symbols():
    return np.array(["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"])


def element_names():
    return np.array(["Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium", "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinum", "Gold", "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", "Francium", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium", "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium", "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium", "Flerovium", "Moscovium", "Livermorium", "Tennessine", "Oganesson"])


# Simulation parameters
dt = 100  # s
time_step = 0
total_time_steps = 100000
total_time = dt*total_time_steps

# Neutron flux parameters
phi = 1e14
fibre_radius = 0.000125

# Material parameters
topes = {}

Ge_concentration = 0
B_concentration = 0
Si_concentration = 1-Ge_concentration-B_concentration
Si28 = Isotope(14, 28, 0.922545/3*Si_concentration)
Si29 = Isotope(14, 29, 0.04672/3*Si_concentration)
Si30 = Isotope(14, 30, 0.030715/3*Si_concentration)
O16 = Isotope(8, 16, 0.99757*(2/3*Si_concentration+2 /
              3*Ge_concentration+3/5*B_concentration))
O17 = Isotope(8, 17, 0.0003838*(2/3*Si_concentration+2 /
              3*Ge_concentration+3/5*B_concentration))
O18 = Isotope(8, 18, 0.002045*(2/3*Si_concentration+2 /
              3*Ge_concentration+3/5*B_concentration))
if Ge_concentration > 0:
    Ge70 = Isotope(32, 70, 0.2052/3*Ge_concentration)
    Ge72 = Isotope(32, 72, 0.2745/3*Ge_concentration)
    Ge73 = Isotope(32, 73, 0.0776/3*Ge_concentration)
    Ge74 = Isotope(32, 74, 0.3652/3*Ge_concentration)
    Ge76 = Isotope(32, 76, 0.0775/3*Ge_concentration)
    # Ge70 = Isotope(32, 70, 0.2052*Ge_concentration)
    # Ge72 = Isotope(32, 72, 0.2745*Ge_concentration)
    # Ge73 = Isotope(32, 73, 0.0776*Ge_concentration)
    # Ge74 = Isotope(32, 74, 0.3652*Ge_concentration)
    # Ge76 = Isotope(32, 76, 0.0775*Ge_concentration)
if B_concentration > 0:
    B10 = Isotope(5, 10, 0.1965*2/5*B_concentration)
    B11 = Isotope(5, 11, 0.8035*2/5*B_concentration)
    # B10 = Isotope(5, 10, 0.1965*B_concentration)
    # B11 = Isotope(5, 11, 0.8035*B_concentration)

if Si_concentration > 0:
    topes[Si28.key] = Si28
    topes[Si29.key] = Si29
    topes[Si30.key] = Si30
topes[O16.key] = O16
topes[O17.key] = O17
topes[O18.key] = O18
if Ge_concentration > 0:
    topes[Ge70.key] = Ge70
    topes[Ge72.key] = Ge72
    topes[Ge73.key] = Ge73
    topes[Ge74.key] = Ge74
    topes[Ge76.key] = Ge76
if B_concentration > 0:
    topes[B10.key] = B10
    topes[B11.key] = B11

# Run the simulation
for i in tqdm(range(total_time_steps)):
    # print(
    #     f"====================================================================================\nTime step: {time_step}")
    # print(topes)
    keys = list(topes.keys())
    for tope in keys:
        topes[tope].N[time_step+1] = topes[tope].N[time_step]

    for tope in keys:
        topes[tope].decay(dt)
        topes[tope].neutron_capture(dt, phi, fibre_radius)
    time_step += 1


# Print the results
for tope in topes:
    topes[tope].N_max = topes[tope].N.max()
    print(f"\n\n{topes[tope].name} has a maximum concentration of {topes[tope].N_max}")
    print(topes[tope])
    print(topes[tope].N)

# Print a list of all tope names
print("All topes present:")
str = ""
for tope in topes:
    str += f"{topes[tope].name}, "
print(str)


relevant_topes = [topes[tope].key for tope in topes if topes[tope].N_max > 1e-19]

# Print a list of all relevant tope names
print("All topes present in concentrations greater than 1e-19:")
str = ""
for tope in relevant_topes:
    str += f"{topes[tope].name}, "
print(str)

# Plot the results
print("\n\n============================\nPlotting...\n")
fig = plt.figure(figsize=(18, 50))
ax = fig.add_subplot(111)
for tope in relevant_topes:
    print(f"Plotting {topes[tope].name}...")
    ax.plot(np.linspace(0, dt*(total_time_steps+1)/3600/24, total_time_steps+1),
            topes[tope].N, label=topes[tope].name)
    ax.annotate(f"{topes[tope].name}", (1.01*total_time/3600/24, 0.85*topes[tope].N[-1]))
ax.set_yscale("log")
plt.ylabel("Proportion of atoms")
plt.xlabel("Time (days)")
plt.title(
    r"Isotope concentration of 3mol.% GeO$_2$ doped silica in neutron flux ($\phi = 10^{14}$ n/cm$^2$/s)")
ax.set_ylim(1e-19, 1)
plt.legend(bbox_to_anchor=(1.08, 0), loc='lower right')
plt.show()


# Calculate the half life of Boron
if B_concentration > 0:
    start_concentration = topes[B10.key].N[0]
    end_concentration = topes[B10.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of B10 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of B10 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")

    start_concentration = topes[B11.key].N[0]
    end_concentration = topes[B11.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of B11 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of B11 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")


# Calculate the half life of Germanium
if Ge_concentration > 0:
    start_concentration = topes[Ge70.key].N[0]
    end_concentration = topes[Ge70.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of Ge70 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of Ge70 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")

    start_concentration = topes[Ge72.key].N[0]
    end_concentration = topes[Ge72.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of Ge72 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of Ge72 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")

    start_concentration = topes[Ge73.key].N[0]
    end_concentration = topes[Ge73.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of Ge73 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of Ge73 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")

    start_concentration = topes[Ge74.key].N[0]
    end_concentration = topes[Ge74.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of Ge74 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of Ge74 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")

    start_concentration = topes[Ge76.key].N[0]
    end_concentration = topes[Ge76.key].N[-1]
    delta_concentration = start_concentration - end_concentration
    fractional_change = delta_concentration/start_concentration

    total_time_days = total_time / 3600 / 24

    print(f"{fractional_change*100}% of Ge76 was burned over {total_time_days} days")

    half_life_B = -np.log(2)/np.log(1-fractional_change)*total_time_days
    print(
        f"The estimated half life of Ge76 in a neutron flux of {phi} n/cm^2/s is {half_life_B} days")


start_population = sum([topes[tope].N[0] for tope in topes])
final_population = sum([topes[tope].N[-1] for tope in topes])
He_population = [topes[tope].N[-1] for tope in topes if topes[tope].name == f"4He"][0]
print(f"Initial population: {start_population}, Final population: {final_population}")
print(f"Helium population: {He_population}")
