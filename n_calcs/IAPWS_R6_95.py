# An implementation of the Internation Association for the Properties of Water and Steam (IAPWS)
# Formulaiton 1995

import numpy as np

Tc = 647.096  # K
rho_c = 322  # kg/m^3
R = 0.46151805  # kJ/(kg K)


def helmholtz_free_energy(T, rho):
    # Calculates the Helmholtz free energy of water using the IAPWS R6-95 model

    # Args: T - temperature in Kelvin
    #       rho - density in kg/m^3
    # Returns: f - Helmholtz free energy in kJ/kg

    f = R*T*(helmholtz_free_energy_phi_r(T, rho) +
             helmholtz_free_energy_phi_0(T, rho))


def helmholtz_free_energy_phi_0(T, rho):
    # Calculates the ideal gas part of the Helmholtz free energy

    # Args: T - temperature in Kelvin
    #       rho - density in kg/m^3
    # Returns: phi_0 - ideal gas part of the Helmholtz free energy in kJ/kg

    delta = rho / rho_c
    tau = Tc / T

    n, gamma = table_1_constants()

    phi_0 = np.log(delta) + n[0] + n[1]*tau + n[2]*np.log(tau)
    for i in range(3, len(n)):
        phi_0 += n[i] * np.log(1-np.exp(-gamma[i]*tau))

    return phi_0


def helmholtz_free_energy_phi_r(T, rho):
    # Calculates the residual part of the Helmholtz free energy

    # Args: T - temperature in Kelvin
    #       rho - density in kg/m^3
    # Returns: A_r - residual part of the Helmholtz free energy in kJ/kg

    tau = Tc / T
    delta = rho / rho_c

    a, b, c, d, t, A, B, C, D, alpha, beta, gamma, epsilon, n = table_2_constants()

    phi_r = 0

    for i in range(0, 6):
        phi_r += n[i] * delta**d[i] * tau**t[i]

    for i in range(7, 50):
        phi_r += n[i] * delta**d[i] * tau**t[i] * np.exp(-delta**c[i])

    for i in range(51, 53):
        phi_r += n[i] * delta**d[i] * tau**t[i] * \
            np.exp(-alpha[i]*(delta-epsilon[i])**2-beta[i]*(tau-gamma[i])**2)

    for i in range(54, 55):
        theta = (1-tau) + A[i]*((delta-B[i])**2) ** (1/(2*beta[i]))
        Delta = theta**2 + B[i]*((delta-1)**2)**a[i]
        psi = np.exp(-C[i]*(delta-1)**2-D[i]*(tau-1)**2)

        phi_r += n[i] * Delta**b[i] * delta * psi

    return phi_r


def helmholtz_free_energy_phi_r_derivative_delta(T, rho):
    # Calculates the derivative of the residual part of the Helmholtz free energy with respect to density

    # Args: T - temperature in Kelvin
    #       rho - density in kg/m^3
    # Returns: phi_r_derivative_delta - derivative of the residual part of the Helmholtz free energy with respect to density in kJ/kg

    tau = Tc / T
    delta = rho / rho_c

    a, b, c, d, t, A, B, C, D, alpha, beta, gamma, epsilon, n = table_2_constants()

    phi_r_derivative_delta = 0

    for i in range(0, 6):
        phi_r_derivative_delta += n[i] * d[i] * delta**(d[i]-1) * tau**t[i]

    for i in range(7, 50):
        phi_r_derivative_delta += n[i] * np.exp(-delta**c[i]) * \
            delta**(d[i]-1) * tau**t[i] * (d[i] - c[i]*delta**c[i])

    for i in range(51, 53):
        phi_r_derivative_delta += n[i] * d[i] * tau**t[i] * np.exp(-alpha[i]*(
            delta-epsilon[i])**2-beta[i]*(tau-gamma[i])**2) * (d[i]/delta - 2*alpha[i]*(delta-epsilon[i]))

    for i in range(54, 55):
        theta = (1-tau) + A[i]*((delta-B[i])**2) ** (1/(2*beta[i]))
        Delta = theta**2 + B[i]*((delta-1)**2)**a[i]
        psi = np.exp(-C[i]*(delta-1)**2-D[i]*(tau-1)**2)

        blergh = (delta-1)*(A[i]*theta*2/beta[i]*((delta-1)**2) **
                            (1/(2*beta[i])-1) + 2*B[i]*a[i]*((delta-1)**2)**(a[i]-1))
        blergh2 = b[i] * Delta**(b[i]-1) * blergh

        blergh3 = -2*C[i]*(delta-1)*psi

        phi_r_derivative_delta += n[i] * \
            (Delta**b[i]*(psi+delta*blergh3)+blergh2*delta*psi)

    return phi_r_derivative_delta


def pressure(T, rho):
    # Calculates the pressure of water using the IAPWS R6-95 model

    # Args: T - temperature in Kelvin
    #       rho - density in kg/m^3
    # Returns: P - pressure in MPa

    delta = rho / rho_c

    P = rho*R*T*(1+delta*helmholtz_free_energy_phi_r_derivative_delta(T, rho))

    return P


def density(P, T):
    # Calculates the density of water using the IAPWS R6-95 model
    # Uses Newton's method to solve for density

    # Args: P - pressure in MPa
    #       T - temperature in Kelvin
    # Returns: rho - density in kg/m^3

    rho_guess = 700  # kg/m^3
    pressure_result = pressure(rho_guess, T)

    delta_rho = 0.01  # kg/m^3

    iterations = 0
    while (np.abs(pressure_result-P) > 0.0000000001 and iterations < 20):
        pressure_result = pressure(rho_guess, T)
        delta_pressure = pressure(
            rho_guess+delta_rho, T) - pressure_result
        rho_guess = rho_guess - (pressure_result-P)/(delta_pressure/delta_rho)
        pressure_result = pressure(rho_guess, T)
        iterations += 1

    return rho_guess


def table_1_constants():
    # Returns the constants for table 1 of the IAPWS R6-95 model

    n = np.array([-8.3204464837497,
                 6.6832105275932,
                 3.00632,
                 0.012436,
                 0.97315,
                 1.27950,
                 0.96956,
                 0.24873])
    gamma = np.array([None,
                      None,
                      None,
                      1.28728967,
                      3.53734222,
                      7.74073708,
                      9.24437796,
                      27.5075105
                      ])
    return n, gamma


def table_2_constants():
    # Returns the constants for table 2 of the IAPWS R6-95 model

    # print(f"rho: {rho_guess}, P: {pressure_result}")
    a = np.concatenate((54*[None], [3.5, 3.5]))

    b = np.concatenate((54*[None], [0.85, 0.95]))

    c = np.concatenate((7*[None],
                       np.ones(15),
                       2*np.ones(20),
                       3*np.ones(4),
                       [4],
                       6*np.ones(4),
                       5*[None]))

    d = np.array([1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2,
                 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6, 3, 3, 3])

    t = np.array([-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1,
                 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50, 0, 1, 4])

    A = np.concatenate((54*[None], [0.32, 0.32]))

    B = np.concatenate((54*[None], [0.2, 0.2]))

    C = np.concatenate((54*[None], [28, 32]))

    D = np.concatenate((54*[None], [700, 800]))

    alpha = np.concatenate((51*[None], [20, 20, 20, None, None]))

    beta = np.concatenate((51*[None], [150, 150, 250, 0.3, 0.3]))

    gamma = np.concatenate((51*[None], [1.21, 1.21, 1.25, None, None]))

    epsilon = np.concatenate((51*[None], [1, 1, 1, None, None]))

    n = np.array([0.12533547935523e-1,
                  0.78957634722828e1,
                  -0.87803203303561e1,
                  0.31802509345418,
                  -0.26145533859358,
                  -0.78199751687981e-2,
                  0.88089493102134e-2,
                  -0.66856572307965,
                  0.20433810950965,
                  -0.66212605039687e-4,
                  -0.19232721156002,
                  -0.25709043003438,
                  0.16074868486251,
                  -0.40092828925807e-1,
                  0.39343422603254e-6,
                  -0.759413770881e-5,
                  0.56250979351888e-3,
                  -0.15608652257135e-4,
                  0.11537996422951e-8,
                  0.36582165144204e-6,
                  -0.13251180074668e-11,
                  -0.62639586912454e-9,
                  -0.10793600908932,
                  0.17611491008752e-1,
                  0.22132295167546,
                  -0.40247669763528,
                  0.58083399985759,
                  0.49969146990806e-2,
                  -0.31358700712549e-1,
                  -0.74315929710341,
                  0.47807329915480,
                  0.20527940895948e-1,
                  -0.13636435110343,
                  0.14180634400617e-1,
                  0.83326504880713e-2,
                  -0.29052536009585e-1,
                  0.38615085574206e-1,
                  -0.20393486513704e-1,
                  -0.16554050063734e-2,
                  0.19955571979541e-2,
                  0.15870308324157e-3,
                  -0.16388568342530e-4,
                  0.43613615723811e-1,
                  0.34994005463765e-1,
                  -0.76788197844621e-1,
                  0.22446277332006e-1,
                  -0.62689710414685e-4,
                  -0.55711118565645e-9,
                  -0.19905718354408,
                  0.31777497330738,
                  -0.11841182425981,
                  -0.31306260323435e2,
                  0.31546140237781e2,
                  -0.25213154341695e4,
                  -0.14874640856724,
                  0.31806110878444])

    return a, b, c, d, t, A, B, C, D, alpha, beta, gamma, epsilon, n


def main():
    # Main function for testing the IAPWS R6-95 model
    print(table_2_constants())


if __name__ == "__main__":
    main()
