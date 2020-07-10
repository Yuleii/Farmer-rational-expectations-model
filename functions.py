""" This implements the solution of a multivariate rational expectations model.

The computations performed in the following code correspond to pages 41 - 51 of
the book by Roger E. A. Farmer, "Macroeconomics of Self-fulfilling
Prophecies", 2nd edition, MIT Press, 1999.

"""
# Import Python Packages
import numpy as np
import matplotlib.pyplot as plt


def solve_farmer_model(alpha, beta, delta, rho, T, inital_productivity_state):
    """Solve a multivariate rational expectations model.

    This function determine the rational expectations solution for consumption
    and capital choices, given realizations of draws of fundamental disturbances
    and the implied time-series of the state of productivity.

    Parameters
    ----------
    alpha: float
        Capital share.
    beta: float
        Discount factor.
    delta: float
        Depreciation rate.
    T: int
        The number of time periods to simulate.
    inital_10_productivity: list
        A list contains 10 inital productivity state values(float).

    Returns
    -------
    s_t: array_like
        Productivity state over time. The shape should be (T, 1).
    c_t: array_like
        Consumption over time. The shape should be (T, 1).
    k_t: array_like
        Capital over time. The shape should be (T+1, 1).
    y_t: array_like
        Output over time. The shape should be (T+1, 1).
    """
    # Getting non-stochastic steady-state values.
    k_bar = ((1 / beta - 1 + delta) / alpha) ** (1 / (alpha - 1))
    c_bar = k_bar ** alpha - delta * k_bar
    # s_bar = 1

    # Describing relationships underlying linear coefficients.
    a1 = beta * alpha * (alpha - 1) * k_bar ** (alpha - 1)
    a2 = beta * alpha * k_bar ** (alpha - 1)
    b1 = 1 - delta + alpha * k_bar ** (alpha - 1)
    b2 = k_bar ** (alpha - 1)
    b3 = -c_bar / k_bar

    # Writing dynamic equilibrium conditions into matrix form.
    # m371 stands for the 1st matrix showing up in equation (3.7), i.e. the
    # matrix multiplying period-t variables on the LHS of (3.7).
    m371 = np.zeros((3, 3))
    m371[0, 0] = -1
    m371[1, 0] = b3
    m371[1, 1] = b1
    m371[1, 2] = b2
    m371[2, 2] = rho

    # m372 stands for the 2nd matrix showing up in equation (3.7)
    m372 = np.zeros((3, 3))
    m372[0, 0] = -1
    m372[0, 1] = a1
    m372[0, 2] = a2
    m372[1, 1] = 1
    m372[2, 2] = 1

    # m38 stands for the inverse of m371
    m38 = np.linalg.inv(m371)

    # Matrix A
    A = np.dot(m38, m372)

    # Calculate Eigenvectors and Eigenvalues of A
    Lam_temp, Q_temp = np.linalg.eig(A)

    # Identify the stable eigenvalue, and therefore can be used to derive a
    # restriction to express a free variable as an equilibrium function of
    # predetermined variables
    select_i_stable = int(np.where(np.absolute(Lam_temp) < 1)[0])
    mark_unstable_1, mark_unstable_2 = np.where(np.absolute(Lam_temp) > 1)[0]

    # Reordering ti make sure the stable eigenvalue and the corresponding
    # eigenvector are in the first positions in the matrix of eigenvalues
    # and the matrix of eigenvectors.
    # Lam = np.diag(
    #     np.array([Q_temp[select_i_stable], Q_temp[mark_unstable_1], Q_temp[mark_unstable_2]]))
    Q = np.zeros((3, 3))
    Q_temp[:, 2]
    Q[:, 0] = Q_temp[:, select_i_stable]
    Q[:, 1] = Q_temp[:, mark_unstable_1]
    Q[:, 2] = Q_temp[:, mark_unstable_2]

    # The inverse of Q, as used in the decomposition A = Q*Lam*inv(Q)
    Q_inv = np.linalg.inv(Q)

    # The key restriction: getting a free variable (consumption) as an
    # equilibrium function of predetermined variables (capital, productivity-state).
    c_coeff_k = -Q_inv[0, 1] / Q_inv[0, 0]
    c_coeff_s = -Q_inv[0, 2] / Q_inv[0, 0]

    # Simulate solution for T periods.
    # Create vector for productivity values.
    s_t = np.zeros((T, 1))

    # First 10 periods given.
    s_t[0:10] = np.array(inital_productivity_state).reshape((-1, 1))

    # Obtain innovations for the productivity process for subsequent periods.
    var_v = 0.007 ** 2
    v_t = np.zeros((T, 1))
    v_t[10::] = np.sqrt(var_v) * np.random.randn(len(v_t) - 10, 1)

    # Obtain productivity levels
    for t in range(10, T):
        s_t[t] = rho * s_t[t - 1] + v_t[t]

    # Generating the time series of variables, which are equilibria for the
    # given shock realizations
    c_t = np.zeros((T, 1))
    k_t = np.zeros((T + 1, 1))

    # Imposing the key equilibrium restriction from rational expectations
    for t in range(0, T):
        c_t[t] = c_coeff_k * k_t[t] + c_coeff_s * s_t[t]
        k_t[t + 1] = b1 * k_t[t] + b2 * s_t[t] + b3 * c_t[t]
    # log linearization of production function y_t = s_t*(k^{alpha}_t)， yields
    # \hat{y}_t = \hat{s}_t + alpha*\hat{k}_t
    y_t = np.transpose(s_t) + alpha * k_t[0:len(k_t) - 1]

    return s_t, c_t, k_t, y_t
