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


def element_symbols():
    return np.array(["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"])


def element_names():
    return np.array(["Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium", "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinum", "Gold", "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", "Francium", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium", "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium", "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium", "Flerovium", "Moscovium", "Livermorium", "Tennessine", "Oganesson"])


class Isotope:
    def __init__(self, Z, A, N):
        self.Z = Z
        self.A = A
        # fraction of atoms in material >>>>>>>>>> NOT NUMBER OF NEUTRONS <<<<<<<<<<
        self.elementSymbol = f"{element_symbols()[Z-1]}"
        self.elementName = f"{element_names()[Z-1]}"
        self.name = f"{A}{self.elementSymbol}"
        self.key = f"{Z}_{A}"
        self.N_start = N
        self.get_parameters()

    def __str__(self):
        str = f"=========================\nElement: {self.elementName}\nIsotope: {self.name}\nZ: {self.Z}\nA: {self.A}"
        if hasattr(self, "half_life"):
            str += f"\nHalf-life: {self.half_life} "
        if hasattr(self, "decay_constant"):
            str += f"\nDecay constant: {self.decay_constant} s^-1"
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
        if hasattr(self, "half_life"):
            total_fraction = 0
            if hasattr(self, "a_fraction"):
                total_fraction += self.a_fraction
            if hasattr(self, "bp_fraction"):
                total_fraction += self.bp_fraction
            if hasattr(self, "bm_fraction"):
                total_fraction += self.bm_fraction
            if hasattr(self, "g_fraction"):
                if total_fraction == 0:
                    total_fraction = self.g_fraction
            if np.abs(1-total_fraction) > 0.01:
                print(
                    f"Error: Decay mode fractions do not add up to 1 ({total_fraction}) - check data for {self.name}")

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
                if row[0] == "0":
                    break
                if line_count == 0:
                    column_half_life = row.index("half_life_sec")
                    column_fraction = row.index("decay_%")
                if line_count == 1:
                    self.half_life = float(row[column_half_life])
                    self.decay_constant = np.log(2)/self.half_life
                    if rad_type == "bm":
                        try:
                            self.bm_fraction = float(row[column_fraction])/100
                        except:
                            self.bm_fraction = 1
                        break
                    elif rad_type == "g":
                        try:
                            self.g_fraction = float(row[column_fraction])/100
                        except:
                            self.g_fraction = 1
                        break
                    elif rad_type == "bp":
                        try:
                            self.bp_fraction = float(row[column_fraction])/100
                        except:
                            self.bp_fraction = 1
                        break
                    elif rad_type == "a":
                        try:
                            self.a_fraction = float(row[column_fraction])/100
                        except:
                            self.a_fraction = 1
                        break
                line_count += 1

    def daughters(self):
        daughters = []
        if hasattr(self, "a_fraction"):
            daughters.append((self.Z-2, self.A-4))
        if hasattr(self, "bp_fraction"):
            daughters.append((self.Z-1, self.A))
        if hasattr(self, "bm_fraction"):
            daughters.append((self.Z+1, self.A))
        if hasattr(self, "n_gamma_xs"):
            daughters.append((self.Z, self.A+1))
        if hasattr(self, "n_alpha_xs"):
            daughters.append((self.Z-2, self.A-3))
            daughters.append((2, 4))
        if hasattr(self, "n_proton_xs"):
            daughters.append((self.Z-1, self.A))

        return daughters


def compute_daughter_isotopes(topes):

    isotope_added = False
    for i in range(len(topes)):
        for daughter in topes[i].daughters():
            not_in_topes = True
            for tope in topes:
                if tope.Z == daughter[0] and tope.A == daughter[1]:
                    not_in_topes = False
            if not_in_topes:
                topes.append(Isotope(daughter[0], daughter[1], 0))
                isotope_added = True

    return topes, isotope_added


def generate_transmutation_matrix(topes, dt, phi):
    # This function generates the matrix for the system of equations

    A = np.identity(len(topes))

    for tope in topes:
        row = topes.index(tope)
        A = decay_and_capture(A, topes, row, dt, phi)
        A = add_to_daughter_population(A, topes, row, dt, phi)

    for tope in topes:
        col = topes.index(tope)
        if np.abs(np.sum(A[:, col])) > 1:
            # print(f"Error: Sum of column {col} is greater than 1")
            A[:, col] /= np.sum(A[:, col])

    # for row in range(len(topes)):
    #     for col in range(len(topes)):
    #         if A[row, col] < 0:
    #             # print(f"Error: Negative value in matrix at {row}, {col} ({A[row, col]})")
    #             # print(topes[row])
    #             # print(f"dt: {dt}")
    #             # print(
    #             #     f"topes[{row}].decay_constant * dt = {topes[row].decay_constant * dt}")
    #             # print(A)
    #         if A[row, col] > 1:
            # print(
            #     f"Error: Value in matrix at {row}, {col} is greater than 1 ({A[row, col]})")
            # print("Parent isotope:")
            # print(topes[col])
            # print("Daughter isotope:")
            # print(topes[row])
            # print(f"dt: {dt}")
            # print(A)

    return A


def decay_and_capture(A, topes, index, dt, phi):
    # This function computes how much an isotope population decreases by

    N = len(topes)
    tope = topes[index]

    # Decay of the index-th isotope
    if hasattr(tope, "decay_constant"):
        A[index, index] -= topes[index].decay_constant * dt

    # Neutron capture by the index-th isotope
    if hasattr(tope, "n_gamma_xs"):
        A[index, index] -= phi * tope.n_gamma_xs * dt
    if hasattr(tope, "n_alpha_xs"):
        A[index, index] -= phi * tope.n_alpha_xs * dt
    if hasattr(tope, "n_proton_xs"):
        A[index, index] -= phi * tope.n_proton_xs * dt

    A[index, index] = max(A[index, index], -1)

    return A


def add_to_daughter_population(A, topes, parent_index, dt, phi):
    # This function adds the population of daughter isotopes to the transmutation matrix

    parent_tope = topes[parent_index]

    for daughter_index in range(len(topes)):
        daughter_tope = topes[daughter_index]
        # Decay of the parent isotope
        if hasattr(parent_tope, "a_fraction"):
            if daughter_tope.Z == parent_tope.Z-2 and daughter_tope.A == parent_tope.A-4:
                A[daughter_index, parent_index] += parent_tope.a_fraction * \
                    parent_tope.decay_constant * dt
            if daughter_tope.Z == 2 and daughter_tope.A == 4:
                A[daughter_index, parent_index] += parent_tope.a_fraction * \
                    parent_tope.decay_constant * dt
        if hasattr(parent_tope, "bp_fraction"):
            if daughter_tope.Z == parent_tope.Z-1 and daughter_tope.A == parent_tope.A:
                A[daughter_index, parent_index] += parent_tope.bp_fraction * \
                    parent_tope.decay_constant * dt
        if hasattr(parent_tope, "bm_fraction"):
            if daughter_tope.Z == parent_tope.Z+1 and daughter_tope.A == parent_tope.A:
                A[daughter_index, parent_index] += parent_tope.bm_fraction * \
                    parent_tope.decay_constant * dt

        # Neutron capture by the parent isotope
        if hasattr(parent_tope, "n_gamma_xs"):
            if topes[daughter_index].Z == parent_tope.Z and topes[daughter_index].A == parent_tope.A+1:
                A[daughter_index, parent_index] += parent_tope.n_gamma_xs * dt * phi
        if hasattr(parent_tope, "n_alpha_xs"):
            if topes[daughter_index].Z == parent_tope.Z-2 and topes[daughter_index].A == parent_tope.A-3:
                A[daughter_index, parent_index] += parent_tope.n_alpha_xs * dt * phi
            if topes[daughter_index].Z == 2 and topes[daughter_index].A == 4:
                A[daughter_index, parent_index] += parent_tope.n_alpha_xs * dt * phi
        if hasattr(parent_tope, "n_proton_xs"):
            if topes[daughter_index].Z == parent_tope.Z-1 and topes[daughter_index].A == parent_tope.A:
                A[daughter_index, parent_index] += parent_tope.n_proton_xs * dt * phi

    return A


def expand_population_matrix(N, number_of_new_isotopes):
    # This function expands the population matrix to accomodate new isotopes

    N_expanded = np.zeros((number_of_new_isotopes, N.shape[1]))

    N_expanded[:N.shape[0], :] = N

    return N_expanded


# Starting Population
Ge_concentration = 0
B_concentration = 0
Si_concentration = 1-Ge_concentration-B_concentration
# Si28 = Isotope(14, 28, 0.922545/3*Si_concentration)
# Si29 = Isotope(14, 29, 0.04672/3*Si_concentration)
# Si30 = Isotope(14, 30, 0.030715/3*Si_concentration)
Si28 = Isotope(14, 28, 0.922545*Si_concentration)
Si29 = Isotope(14, 29, 0.04672*Si_concentration)
Si30 = Isotope(14, 30, 0.030715*Si_concentration)
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
if B_concentration > 0:
    B10 = Isotope(5, 10, 0.1965*2/5*B_concentration)
    B11 = Isotope(5, 11, 0.8035*2/5*B_concentration)


topes = [Si28, Si29, Si30, O16, O17, O18]
# topes = [Si28, Si29, Si30]
if Ge_concentration > 0:
    topes += [Ge70, Ge72, Ge73, Ge74, Ge76]
if B_concentration > 0:
    topes += [B10, B11]

# Simulation Parameters
total_time_steps = 100000
dt = 100  # s

# Neutron flux
phi = 1e14  # n/cm^2/s

# Population matrix
N = np.zeros((len(topes), total_time_steps))
for tope in range(len(topes)):
    N[tope, 0] = topes[tope].N_start

# Transmutation matrix
A = generate_transmutation_matrix(topes, dt, phi)

# Run simulation
for time_step in tqdm(range(1, total_time_steps)):
    # print("Time step: ", time_step)

    # Calculate daughter isotopes
    topes, isotope_added = compute_daughter_isotopes(topes)

    # Regenerate matrix if new isotopes were added
    if isotope_added:
        A = generate_transmutation_matrix(topes, dt, phi)
        N = expand_population_matrix(N, len(topes))

    # Calculate new population
    N[:, time_step] = np.matmul(A, N[:, time_step-1])

# Calculate the maximum population of each isotope and print the isotope details.
for tope in range(len(topes)):
    topes[tope].N_max = np.max(N[tope, :])
    print(topes[tope])

# Generate a list of isotopes with a significant population
relevant_tope_indices = [i for i in range(len(topes)) if topes[i].N_max > 1e-19]

# Plotting
dt_days = dt/3600/24
time = np.linspace(0, dt_days*total_time_steps, total_time_steps)


fig = plt.figure(figsize=(18, 50))
ax = fig.add_subplot(111)
for i in relevant_tope_indices:
    ax.plot(time, N[i, :], label=topes[i].name)
    ax.annotate(f"{topes[i].name}", (1.01*time[-1], 0.85*N[i, -1]))
ax.set_yscale("log")
ax.set_xscale("log")
plt.ylabel("Proportion of atoms")
plt.xlabel("Time (days)")
plt.grid()
plt.title(
    r"Isotope concentration of pure silica glass in neutron flux ($\phi = 10^{14}$ n/cm$^2$/s)")
ax.set_ylim(1e-19, 1)
plt.legend(bbox_to_anchor=(1.08, 0), loc='lower right')
plt.show()

print(f"\nStarting population: {sum(N[:, 0])}, Final population: {sum(N[:, -1])}")
He_index = [i for i in range(len(topes)) if topes[i].name == "4He"][0]
print(f"Helium-4 population: {N[He_index, -1]}")


fig = plt.figure(figsize=(18, 20))
plt.imshow(np.log10(A))
plt.colorbar()
plt.show()
