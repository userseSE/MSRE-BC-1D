from scipy.sparse.linalg import spsolve
import numpy as np
from scipy.sparse import csc_matrix

from ode_solver import ode_solver


def neutronics(y_n, rho, step, params):
    # Extract relevant parameters from the params dictionary
    v_core = params['v_core']
    inlet_mode = params['inlet_mode']
    N = params['N']
    dt = params['dt']
    dz = params['dz']
    V = params['V']
    D = params['D']
    sigma_a = params['sigma_a']
    nu_sigma_f = params['nu_sigma_f']

    beta = params.get('beta_np', np.asarray(params['beta'], dtype=float))          # array-like, len=6
    lambda_i = params.get('lambda_i_np', np.asarray(params['lambda_i'], dtype=float))  # array-like, len=6
    Beta = float(beta.sum())

    # --- NEW (minimal): drift settings ---
    v_core = params.get('v_core', 0.2)              # m/s
    inlet_mode = params.get('inlet_mode', 'copy')   # 'fresh' or 'copy'
    if inlet_mode not in ('fresh', 'copy'):
        raise ValueError("params['inlet_mode'] must be 'fresh' or 'copy'")

    A = params.get('A_neutronics')
    B = params.get('B_neutronics')
    if A is None or B is None:
        # Discretize the spatial domain and set up the Crank-Nicolson matrices
        main_diag = -2 * np.ones(N)
        off_diag = np.ones(N - 1)
        D2 = (np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)) / dz**2

        # Dirichlet boundary conditions for zero flux at boundaries
        D2[0, :] = 0
        D2[-1, :] = 0
        D2[0, 0] = 1
        D2[-1, -1] = 1

        D2_sparse = csc_matrix(D2)
        I = np.eye(N)

        A = csc_matrix(I - 0.5 * dt * V * D * D2_sparse)
        B = csc_matrix(I + 0.5 * dt * V * D * D2_sparse)

    # NOTE: your original code uses Keff = 1/(1 - rho) indirectly via (1+rho).
    # Keep your original convention to avoid breaking your model.
    # (If you later want to fix consistency, we can do it separately.)

    def pde_to_ode_neutronics(t, y):
        phi = y[:N]

        # delayed source term sum_i lambda_i * Ci
        C = y[N:].reshape(6, N)
        lambda_ci = (lambda_i[:, None] * C).sum(axis=0)

        # Crank-Nicolson update for phi
        rhs_phi = B @ phi + dt * V * ((-sigma_a + ((1 - Beta) * nu_sigma_f * (1 + rho))) * phi + lambda_ci)
        phi_new = spsolve(A, rhs_phi)
        dphi_dt = (phi_new - phi) / dt

        # ---------- NEW: precursor update with upwind drift ----------
        # Eq: dCi/dt + d/dz(v Ci) = beta_i * nuSigF*(1+rho)*phi - lambda_i*Ci
        # => dCi/dt = production - lambda_i*Ci - d/dz(v Ci)
        # upwind derivative for v>0: d/dz(vC) ~ v*(C[j]-C[j-1])/dz
        # (If v<0, you should switch to downwind; but MSRE core flow is typically positive z here.)
        v = float(v_core)
        C_im1 = np.empty_like(C)
        C_im1[:, 1:] = C[:, :-1]
        if inlet_mode == 'fresh':
            C_im1[:, 0] = 0.0
        else:  # 'copy'
            C_im1[:, 0] = C[:, 0]

        adv = v * (C - C_im1) / dz
        production = (beta[:, None] * (nu_sigma_f * (1 + rho)) * phi)
        dci_dt = production - (lambda_i[:, None] * C) - adv

        # return concatenated derivative vector
        return np.concatenate([dphi_dt, dci_dt.reshape(-1)])

    # Initial condition vector
    if step == 0:
        phi_0 = params['phi_0']
        c0 = params['c0']
        y0 = np.concatenate([
            phi_0,
            c0[:N], c0[N:2 * N], c0[2 * N:3 * N],
            c0[3 * N:4 * N], c0[4 * N:5 * N], c0[5 * N:]
        ])
    else:
        y0 = y_n

    # Solve the system of ODEs
    bc = []
    solution_y_n = ode_solver(y0, bc, pde_to_ode_neutronics, params)

    # Extract the solution
    phi = solution_y_n[:N, -1].T
    q_prime = phi
    y_n = solution_y_n

    return y_n, q_prime
