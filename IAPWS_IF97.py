# This seems to be a very complicated method. The 1995 Formulation seems much simpler...

import numpy as np

R = 0.461526  # kJ/(kg K) - specific gas constant of water


def specific_volume(P, T):
    # Determine the specific volume of water at a given pressure and temperature using the
    # IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: nu - specific volume in m^3/kg

    region = density_region(P, T)

    if region == -1:
        return -1

    if region == 1:
        return region1_specific_volume(P, T)
    if region == 2:
        return region2_specific_volume(P, T)
    if region == 3:
        return region3_specific_volume(P, T)
    if region == 4:
        return -1e16
        # return region4_specific_volume(P, T)
    if region == 5:
        return -1e16
        # return region5_specific_volume(P, T)


def density(P, T):
    # Determine the density of water at a given pressure and temperature using the
    # IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: rho - density in kg/m^3

    nu = specific_volume(P, T)

    rho = 1/nu
    return rho


def region1_specific_volume(P, T):
    # Calculate the specific volume of water in region 1 using the IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: nu - specific volume in m^3/kg

    Pstar = 16.53  # MPa

    pi = P/Pstar

    gamma_pi = region1_gamma_pi(P, T)

    return pi*gamma_pi*R*T/P / 1000


def region1_gamma_pi(P, T):
    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free energy
    # for water in region 1 using the IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: gamma_pi - partial derivative with respect to pi of the dimensionless Gibbs free energy

    n, I, J = table_A2_constants()

    Pstar = 16.53  # MPa
    Tstar = 1386  # K

    pi = P/Pstar
    tau = Tstar/T

    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free energy
    gamma_pi = 0
    for i in range(len(n)):
        gamma_pi += -n[i] * I[i] * (7.1-pi)**(I[i]-1) * (tau-1.222)**J[i]

    return gamma_pi


def region2_specific_volume(P, T):
    # Calculate the specific volume of water in region 2 using the IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: nu - specific volume in m^3/kg

    # Calculate the partial derivative with respect to pi of the deimensionless Gibbs free
    # energy for the ideal gas part
    gamma0_pi = region2_gamma0_pi(P, T)

    Pstar = 1  # MPa
    pi = P/Pstar

    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free
    # energy for the residual part
    gammar_pi = region2_gammar_pi(P, T)

    # Calculate the specific volume
    nu = pi*(gamma0_pi+gammar_pi)*R*T/P / 1000  # m*3/kg
    return nu


def region2_gamma0_pi(P, T):
    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free energy
    # for water in region 2 using the IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: gamma0_pi - partial derivative with respect to pi of the dimensionless Gibbs free energy

    Pstar = 1  # MPa

    pi = P/Pstar

    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free energy
    gamma0_pi = 1/pi
    return gamma0_pi


def region2_gammar_pi(P, T):
    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free energy
    # for water in region 2 using the IAPWS-IF97 equations of state

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: gammar_pi - partial derivative with respect to pi of the dimensionless Gibbs free energy

    Pstar = 1  # MPa
    Tstar = 540  # K

    pi = P/Pstar
    tau = Tstar/T

    nr, Ir, Jr = table_A5_constants()

    # Calculate the partial derivative with respect to pi of the dimensionless Gibbs free
    # energy for the residual part
    gammar_pi = 0
    for i in range(len(nr)):
        gammar_pi += nr[i] * Ir[i] * pi**(Ir[i]-1) * (tau-0.5)**Jr[i]

    return gammar_pi


def region3_specific_volume(P, T):
    return region3_density(P, T)**-1


def region3_density(P, T):
    # Calculate the density of water in region 3 using the IAPWS-IF97 equations of state
    # Since the equation of state takes the density as an input, the density
    # for a given P and T is calculated by guessing and using Newton's method to find the
    # density

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: rho - density in kg/m^3

    if T < 647.096:  # Critical point temperature
        if P > saturation_pressure(T):
            rho_guess = 700  # kg/m^3
        else:
            rho_guess = 0.001
    else:
        if P > auxilliary_pressure(T):
            rho_guess = 700
        else:
            rho_guess = 0.001

    pressure_result = region3_pressure(rho_guess, T)

    delta_rho = 0.01  # kg/m^3

    iterations = 0
    while (np.abs(pressure_result-P) > 0.00000001 and iterations < 20):
        # print(f"rho: {rho_guess}, P: {pressure_result}")

        pressure_result = region3_pressure(rho_guess, T)
        delta_pressure = region3_pressure(
            rho_guess+delta_rho, T) - pressure_result
        rho_guess = rho_guess - (pressure_result-P)/(delta_pressure/delta_rho)
        pressure_result = region3_pressure(rho_guess, T)
        iterations += 1

    # print(f"rho: {rho_guess}, P: {pressure_result}")
    return rho_guess


def region3_pressure(rho, T):
    # Calculate the pressure of water in region 3 using the IAPWS-IF97 equations of state

    # Args: rho - density in kg/m^3
    #       T - temperature in Kelvin
    # Returns: P - pressure in MPa

    rhostar = 322  # kg/m^3

    delta = rho/rhostar

    phi_delta = region3_phi_delta(rho, T)

    P = rho*R*T * delta * phi_delta / 1000  # MPa
    return P


def region3_phi_delta(rho, T):
    # Calculate the dimensionless Helmholtz free energy for water in region 3 using the
    # IAPWS-IF97 equations of state

    # Args: rho - dimensionless density
    #       T - temperature in Kelvin
    # Returns: phi - dimensionless Helmholtz free energy

    n, I, J = table_A9_constants()

    Tstar = 647.096  # K
    rhostar = 322  # kg/m^3

    tau = Tstar/T
    delta = rho/rhostar

    phi_delta = n[0] / delta
    for i in range(1, len(n)):
        phi_delta += n[i] * I[i] * (delta)**(I[i]-1) * (tau)**J[i]

    return phi_delta


def density_region(P, T):
    # Determine the region defined by IAPWS-IF97 for a given pressure and temperature.
    # This determines which equations of state to use to determine density.

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin

    # Returns: region - the region defined by IAPWS-IF97 for the given pressure and temperature

    # Determine if the given pressure and temperature are within the valid range
    if P > 100:
        return -1
    if P < 0:
        return -1
    if T > 1073.15 and P > 50:
        return -1
    if T > 2273.15:
        return -1
    if T < 0:
        return -1

    # A small range of pressure is used to determine if the pressure is on the saturation line
    deltaP = 0.000001  # MPa
    if T > 1073.15:
        return 5
    if T > 623.15:
        if P < auxilliary_pressure(T):
            return 2
        else:
            return 3
    else:
        if np.abs(P - saturation_pressure(T)) < deltaP:
            return 4
        if P < saturation_pressure(T):
            return 2
        else:
            return 1


def saturation_pressure(T):
    # Calculate the saturation pressure of water at a given temperature using the IAPWS-IF97
    # equations of state

    # Args: T - temperature in Kelvin
    # Returns: P - saturation pressure in MPa

    # Calculate the dimensionless temperature
    theta = dimensionless_temperature(T)

    Pstar = 1  # MPa

    n = table_A11_constants()

    A = theta**2 + n[0]*theta + n[1]
    B = n[2]*theta**2 + n[3]*theta + n[4]
    C = n[5]*theta**2 + n[6]*theta + n[7]

    P = Pstar * (2*C/(-B+(B**2-4*A*C)**0.5))**4
    return P


def dimensionless_temperature(T):
    # Calculate the dimensionless temperature for use in the IAPWS-IF97 equations of state

    # Args: T - temperature in Kelvin
    # Returns: theta - dimensionless temperature

    Tstar = 1

    n = table_A11_constants()

    theta = T/Tstar + n[8]/((T/Tstar)-n[9])
    return theta


def auxilliary_pressure(T):
    # Calculate the auxilliary pressure for use in the IAPWS-IF97 equations of state.
    # This determines the boundary between regions 2 and 3

    # Args: T - temperature in Kelvin
    # Returns: P - auxilliary pressure in MPa

    Tstar = 1  # K
    theta = T/Tstar  # Calculate the dimensionless temperature

    # Constants from table A1 in IAPWS-IF97
    n1 = 0.34805185628969e3
    n2 = -0.11671859879975e1
    n3 = 0.10192970039326e-2

    Pstar = 1  # MPa

    P = Pstar * (n1 + n2*theta + n3*theta**2)
    return P


def table_A2_constants():
    # Constants from table A2 in IAPWS-IF97
    n = np.array([0.14632971213167,
                  -0.84548187169114,
                  -3.7563603672040,
                  3.3855169168385,
                  -0.95791963387872,
                  0.15772038513228,
                  -0.16616417199501e-1,
                  0.81214629983568e-3,
                  0.28319080123804e-3,
                  -0.60706301565874e-3,
                  -0.18990068218419e-1,
                  -0.32529748770505e-1,
                  -0.21841717175414e-1,
                  -0.52838357969930e-4,
                  -0.47184321073267e-3,
                  -0.30001780793026e-3,
                  0.47661393906987e-4,
                  -0.44141845330846e-5,
                  -0.72694996297594e-15,
                  -0.31679644845054e-4,
                  -0.28270797985312e-5,
                  -0.85205128120103e-9,
                  -0.22425281908000e-5,
                  -0.65171222895601e-6,
                  -0.14341729937924e-12,
                  -0.40516996860117e-6,
                  -0.12734301741641e-8,
                  -0.17424871230634e-9,
                  -0.68762131295531e-18,
                  0.14478307828521e-19,
                  0.26335781662795e-22,
                  -0.11947622640071e-22,
                  0.18228094581404e-23,
                  -0.93537087292458e-25])
    I = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                  2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32])
    J = np.array([-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3,
                  17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41])

    return n, I, J


def table_A5_constants():
    # Constants from table A5 in IAPWS-IF97
    nr = np.array([-0.17731742473213e-2,
                   -0.17834862292358e-1,
                   -0.45996013696365e-1,
                   -0.57581259083432e-1,
                   -0.50325278727930e-1,
                   -0.33032641670203e-4,
                   -0.18948987516315e-3,
                   -0.39392777243355e-2,
                   -0.43797295650573e-1,
                   -0.26674547914087e-4,
                   0.20481737692309e-7,
                   0.43870667284435e-6,
                   -0.32277677238570e-4,
                   -0.15033924542148e-2,
                   -0.40668253562649e-1,
                   -0.78847309559367e-9,
                   0.12790717852285e-7,
                   0.48225372718507e-6,
                   0.22922076337661e-5,
                   -0.16714766451061e-10,
                   -0.21171472321355e-2,
                   -0.23895741934104e2,
                   -0.59059564324270e-17,
                   -0.12621808899101e-5,
                   -0.38946842435739e-1,
                   0.11253211360459e-10,
                   -0.82311340897998e1,
                   0.19809712802088e-7,
                   0.10406965210174e-18,
                   -0.10234747095929e-12,
                   -0.10018179379511e-8,
                   -0.80882908646985e-10,
                   0.10693031879409,
                   -0.33662250574171,
                   0.89185845355421e-24,
                   0.30629316876232e-12,
                   -0.42002467698208e-5,
                   -0.59056029685639e-25,
                   0.37826947613457e-5,
                   -0.12768608934681e-14,
                   0.73087610595061e-28,
                   0.55414715350778e-16,
                   -0.94369707241210e-6])
    Ir = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7,
                   7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24])
    Jr = np.array([0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0,
                   11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58])

    return nr, Ir, Jr


def table_A9_constants():
    # Constants from table A9 in IAPWS-IF97
    n = np.array([0.10658070028513e1,
                  -0.15732845290239e2,
                  0.20944396974307e2,
                  -0.76867707878716e1,
                  0.26185947787954e1,
                  -0.28080781148620e1,
                  0.12053369696517e1,
                  -0.84566812812502e-2,
                  -0.12654315477714e1,
                  -0.11524407806681e1,
                  0.88521043984318,
                  -0.64207765181607,
                  0.38493460186671,
                  -0.85214708824206,
                  0.48972281541877e1,
                  -0.30502617256965e1,
                  0.39420536879154e-1,
                  0.12558408424308,
                  -0.27999329698710,
                  0.13899799569460e1,
                  -0.20189915023570e1,
                  -0.82147637173963e-2,
                  -0.47596035734923,
                  0.43984074473500e-1,
                  -0.44476435428739,
                  0.90572070719733,
                  0.70522450087967,
                  0.10770512626332,
                  -0.32913623258954,
                  -0.50871062041158,
                  -0.22175400873096e-1,
                  0.94260751665092e-1,
                  0.16436278447961,
                  -0.13503372241348e-1,
                  -0.14834345352472e-1,
                  0.57922953628084e-3,
                  0.32308904703711e-2,
                  0.80964802996215e-4,
                  -0.16557679795037e-3,
                  -0.44923899061815e-4
                  ])
    I = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3,
                 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11])
    J = np.array([0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26,
                 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26])

    return n, I, J


def table_A11_constants():
    # Constants from table A11 in IAPWS-IF97
    n = np.array([0.11670521452767e4,
                 -0.72421316703206e6,
                 -0.17073846940092e2,
                 0.12020824702470e5,
                 -0.32325550322333e7,
                 0.14915108613530e2,
                 -0.48232657361591e4,
                 0.40511340542057e6,
                 -0.23855557567849,
                 0.65017534844798e3])
    return n
