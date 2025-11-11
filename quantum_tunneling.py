import numpy as np
import math
import matplotlib.pyplot as plt
from numba import njit

hbar = 1.054_571_817e-34     # J*s
m_e  = 9.109_383_7015e-31    # kg, electron mass
eV   = 1.602_176_634e-19     # J

V0_eV      = 1.0            # barrier height (eV)
x_min_a    = -40.0          # left boundary of computational domain (unit: a)
x_max_a    =  80.0          # right boundary of computational domain (unit: a)
Nx         = 6000           # number of grid points (more is more accurate but slower)
dt_dimless = 1.0e-4         # dimensionless time step Δt/t0
absorb_frac= 0.12           # absorbing layer width as fraction of domain (0~0.4 is reasonable)
absorb_pow = 6              # power for absorbing mask; larger means steeper decay
save_fig   = True           # whether to save figures to files

E0_target_eV = 1.0/100.0        # = 0.01 eV
a_m  = math.sqrt(hbar**2 / (2*m_e*E0_target_eV*eV))
a_nm = a_m * 1e9
E0_J  = hbar**2/(2*m_e*a_m**2)
E0_eV = E0_J / eV
t0_s  = hbar / E0_J
V0_tilde = V0_eV / E0_eV def make_grid_and_potential(x_min_a, x_max_a, Nx, V0_tilde):
    x = np.linspace(x_min_a, x_max_a, Nx)         # dimensionless coordinates (unit: a)
    dx = x[1] - x[0]
    V = np.zeros_like(x)
    # rectangular barrier: V0 for 0 <= x <= 1, zero elsewhere
    mask_barrier = (x >= 0.0) & (x < 1.0)
    V[mask_barrier] = V0_tilde
    return x, dx, V, mask_barrier

x, dx, V, mask_barrier = make_grid_and_potential(x_min_a, x_max_a, Nx, V0_tilde) def build_CN_tridiagonals(V, dx, dt):
    """
    Build tridiagonal system for interior DOF only (size N-2):
      H_diag =  2/dx^2 + V_i
      H_off  = -1/dx^2
    Returns (A_lower, A_diag, A_upper), (B_lower, B_diag, B_upper)
    """
    N = V.size
    if N < 3:
        raise ValueError("Grid too small.")
    dx2 = dx*dx
    Vi = V[1:-1]                              # interior point potentials
    H_diag = (2.0/dx2) + Vi
    H_off  = -1.0/dx2

    dtfac = 0.5j*dt
    # A = I + i dt/2 H
    A_diag = 1.0 + dtfac * H_diag
    A_off  = dtfac * H_off
    # B = I - i dt/2 H
    B_diag = 1.0 - dtfac * H_diag
    B_off  = -dtfac * H_off

    n = N-2
    A_lower = np.full(n-1, A_off, dtype=np.complex128)
    A_upper = np.full(n-1, A_off, dtype=np.complex128)
    B_lower = np.full(n-1, B_off, dtype=np.complex128)
    B_upper = np.full(n-1, B_off, dtype=np.complex128)

    return (A_lower, A_diag.astype(np.complex128), A_upper), \
           (B_lower, B_diag.astype(np.complex128), B_upper) def tridiag_mv(lower, diag, upper, v):
    """Tridiagonal matrix-vector multiplication (consistent sizes: diag=n, lower/upper=n-1, v=n)."""
    n = diag.size
    out = diag * v
    out[1:] += lower * v[:-1]
    out[:-1]+= upper * v[1:]
    return out

def tridiag_solve(lower, diag, upper, rhs):
    """Thomas algorithm to solve tridiagonal A x = rhs. Inputs are copied, original arrays unchanged."""
    n = diag.size
    c = upper.copy()
    d = diag.copy()
    b = rhs.copy()

    # forward elimination
    for i in range(1, n):
        w = lower[i-1] / d[i-1]
        d[i] -= w * c[i-1]
        b[i] -= w * b[i-1]

    # back substitution
    x = np.empty_like(rhs)
    x[-1] = b[-1] / d[-1]
    for i in range(n-2, -1, -1):
        x[i] = (b[i] - c[i]*x[i+1]) / d[i]
    return x def gaussian_packet(x, x0, sigma, k0):
    """Normalized Gaussian wave packet (dimensionless): exp(-(x-x0)^2/(4*sigma^2)) * exp(i k0 (x-x0))"""
    psi = np.exp(-0.5*((x-x0)/sigma)**2) * np.exp(1j*k0*(x-x0))
    # normalization
    dx = x[1]-x[0]
    norm = np.sqrt(np.sum(np.abs(psi)**2)*dx)
    return psi / norm

def make_absorb_mask(x, frac=0.12, pow_k=6):
    """
    Create absorbing layers on left and right edges, each taking 'frac' of total length,
    using a smooth cosine-power window. Center region=1, boundaries=0. Larger pow_k means steeper decay.
    """
    N = x.size
    mask = np.ones(N, dtype=np.float64)
    width = int(max(2, frac*N))
    # left absorbing layer
    for i in range(width):
        s = (i+1)/width  # 0->1
        mask[i] = math.sin(0.5*math.pi*s)**pow_k
    # right absorbing layer
    for i in range(width):
        s = (i+1)/width
        mask[-(i+1)] = math.sin(0.5*math.pi*s)**pow_k
    return mask @njit(cache=True)
def thomas_prefactor(lower, diag, upper):
    """
    One-time LU decomposition of tridiagonal matrix A (Thomas forward elimination coefficients).
    Returns:
      w  : multipliers from L factor (length n-1)
      d  : diagonal of U (length n)
      c  : superdiagonal of U (length n-1)
    """
    n = diag.size
    c = upper.copy()
    d = diag.copy()
    w = np.empty(n-1, dtype=diag.dtype)
    for i in range(1, n):
        wi = lower[i-1] / d[i-1]
        w[i-1] = wi
        d[i]  -= wi * c[i-1]
    return w, d, c

@njit(cache=True)
def thomas_solve_fact(w, d, c, rhs):
    """
    Solve A x = rhs using pre-factored coefficients (w,d,c).
    Only forward and backward substitution on RHS, no repeated matrix computation.
    """
    n = d.size
    b = rhs.copy()
    for i in range(1, n):
        b[i] -= w[i-1] * b[i-1]
    x = np.empty_like(rhs)
    x[-1] = b[-1] / d[-1]
    for i in range(n-2, -1, -1):
        x[i] = (b[i] - c[i] * x[i+1]) / d[i]
    return x def run_one_energy(E_eV, x, dx, V, dt_dimless,
                   x0=-12.0, sigma=2.5, store_stride=200,
                   absorb_frac=0.12, absorb_pow=6, t_safety=1.5,
                   stop_tol=1e-4, min_steps=2000, _converge_window=5):
    """
    Crank-Nicolson time evolution for a single energy (with pre-factorization + convergence early stopping).

    Parameters:
      stop_tol : float  convergence threshold (if T and norm both change less than this value
                        over _converge_window records, stop early)
      min_steps: int    minimum number of steps before allowing early stopping
      _converge_window: window size for comparing convergence between consecutive records

    Returns dict with:
      {
        'E_eV': E_eV, 'k0': k0, 'T': T_last, 'R': R_last,
        'T_series': T_series, 'R_series': R_series,
        'norm_series': norms, 'times': times,
        'psi_final': psi, 'x': x
      }
    """
    # 1) Dimensionless energy and wave number
    E_tilde = E_eV / E0_eV
    k0 = math.sqrt(max(E_tilde, 1e-14))  # prevent numerical issues for very small energies

    # 2) Initial state (normalized Gaussian wave packet)
    psi = gaussian_packet(x, x0=x0, sigma=sigma, k0=k0).astype(np.complex128)

    # 3) CN tridiagonal matrices (interior degrees of freedom)
    (Al, Ad, Au), (Bl, Bd, Bu) = build_CN_tridiagonals(V, dx, dt_dimless)

    # 4) Pre-factorize A (use Numba version if available; otherwise fall back to regular solving)
    use_prefactor = False
    try:
        # Only enable if the functions are defined globally
        if 'thomas_prefactor' in globals() and 'thomas_solve_fact' in globals():
            wA, dA, cA = thomas_prefactor(Al, Ad, Au)  # one-time decomposition
            use_prefactor = True
        else:
            wA = dA = cA = None
    except Exception:
        # auto fallback if error occurs
        use_prefactor = False
        wA = dA = cA = None

    # 5) Absorbing mask
    absorb = make_absorb_mask(x, frac=absorb_frac, pow_k=absorb_pow)

    # 6) Estimate maximum evolution steps: from x0 to right side of barrier + buffer, group velocity vg=2k0
    dist   = (x.max() - 1.2) - x0
    vg     = max(2.0*k0, 1e-8)
    t_total= t_safety * dist / vg
    nsteps = int(max(1, t_total / dt_dimless))

    # 7) Monitoring quantities
    rights, lefts, norms, times = [], [], [], []

    # 8) Statistics region masks
    right_mask = (x > 1.05)
    left_mask  = (x < -0.05)

    # 9) Time evolution main loop
    psi_in = psi[1:-1]        # store only interior degrees of freedom
    last_store_step = -1

    for n in range(nsteps):
        # right-hand side vector
        rhs = tridiag_mv(Bl, Bd, Bu, psi_in)

        # solve A * psi_in^{n+1} = rhs
        if use_prefactor:
            psi_in = thomas_solve_fact(wA, dA, cA, rhs)
        else:
            psi_in = tridiag_solve(Al, Ad, Au, rhs)

        # write back + Dirichlet boundary condition
        psi[1:-1] = psi_in
        psi[0]    = 0.0 + 0.0j
        psi[-1]   = 0.0 + 0.0j

        # absorption (apply each step)
        psi *= absorb

        # record (every stride or at last step)
        if (n % store_stride) == 0 or n == nsteps - 1:
            prob = np.abs(psi)**2
            norms.append(prob.sum()*dx)
            rights.append(prob[right_mask].sum()*dx)
            lefts.append(prob[left_mask].sum()*dx)
            times.append((n+1)*dt_dimless)
            last_store_step = n

            # 10) Convergence early stopping: both T and norm change are small, and min_steps exceeded
            if len(rights) > _converge_window and (n+1) > min_steps:
                dT   = abs(rights[-1] - rights[-1 - _converge_window])
                dNrm = abs(norms[-1]  - norms[-1  - _converge_window])
                if dT < stop_tol and dNrm < stop_tol:
                    break

    # ensure at least one record; if loop didn't trigger record, force one record
    if last_store_step < 0:
        prob = np.abs(psi)**2
        norms.append(prob.sum()*dx)
        rights.append(prob[right_mask].sum()*dx)
        lefts.append(prob[left_mask].sum()*dx)
        times.append(dt_dimless)

    # 11) Summarize results
    T_last = float(rights[-1])
    R_last = float(lefts[-1])
    return {
        'E_eV': E_eV, 'k0': k0, 'T': T_last, 'R': R_last,
        'T_series': np.array(rights), 'R_series': np.array(lefts),
        'norm_series': np.array(norms), 'times': np.array(times),
        'psi_final': psi, 'x': x
    } # construct grid and potential
x, dx, V, mask_barrier = make_grid_and_potential(x_min_a, x_max_a, Nx, V0_tilde)

# test a single energy point (can be adjusted)
E_test = 0.8  # eV
res = run_one_energy(E_test, x, dx, V, dt_dimless,
                     x0=-14.0, sigma=3.0,
                     absorb_frac=absorb_frac, absorb_pow=absorb_pow, t_safety=1.6)

print(f"E = {res['E_eV']:.3f} eV,  k0(dimless)={res['k0']:.4f},  "
      f"T≈{res['T']:.4f}, R≈{res['R']:.4f}, T+R≈{res['T']+res['R']:.4f}")

# plot monitoring quantities
plt.figure(figsize=(6,4))
plt.plot(res['times'], res['T_series'], label='T(t) right side probability')
plt.plot(res['times'], res['R_series'], label='R(t) left side probability')
plt.plot(res['times'], res['norm_series'], label='norm ∫|ψ|²dx')
plt.xlabel('t / t0')
plt.ylabel('value')
plt.legend()
plt.tight_layout()
if save_fig: plt.savefig('monitor_single_energy.png', dpi=160)
plt.show()

# plot final spatial distribution
plt.figure(figsize=(7,4))
plt.plot(x, np.abs(res['psi_final'])**2, label='|ψ(x,t_end)|²')
# mark the barrier region (0~1)
plt.axvspan(0.0, 1.0, alpha=0.2, label='barrier region [0,1]')
plt.xlabel('x / a')
plt.ylabel('|ψ|²')
plt.legend()
plt.tight_layout()
if save_fig: plt.savefig('psi_final_single_energy.png', dpi=160)
plt.show() def k_from_E_eV(E_eV, m=m_e):
    """Wave number k (1/m) in free region, input in eV -> J"""
    E_J = E_eV * eV
    return np.sqrt(2*m*E_J)/hbar

def T_analytic_rect_barrier(E_eV, V0_eV=V0_eV, a_m=a_m, m=m_e):
    E = np.array(E_eV, dtype=float)
    V0 = float(V0_eV)
    T = np.empty_like(E, dtype=float)

    above = (E > V0)
    below = ~above

    # E > V0:   sin
    if np.any(above):
        k2 = np.sqrt(2*m*((E[above]-V0)*eV))/hbar
        num = V0**2
        den = 4.0*E[above]*(E[above]-V0)
        T[above] = 1.0 / (1.0 + (num/den) * np.sin(k2*a_m)**2)

    # E < V0:   sinh
    if np.any(below):
        kappa = np.sqrt(2*m*((V0-E[below])*eV))/hbar
        num = V0**2
        den = 4.0*E[below]*(V0-E[below])
        T[below] = 1.0 / (1.0 + (num/den) * np.sinh(kappa*a_m)**2)

    return T  def sweep_energies(E_list_eV, x, dx, V, dt_dimless,
                   x0=-14.0, sigma=3.0,
                   absorb_frac=0.12, absorb_pow=6, t_safety=1.6):
    Tnum = []
    Rnum = []
    for E in E_list_eV:
        res = run_one_energy(E, x, dx, V, dt_dimless,
                             x0=x0, sigma=sigma,
                             absorb_frac=absorb_frac, absorb_pow=absorb_pow, t_safety=t_safety)
        Tnum.append(res['T'])
        Rnum.append(res['R'])
    return np.array(Tnum), np.array(Rnum)

# energy sampling
E_min, E_max, E_num = 0.05, 2.0, 50
E_list = np.linspace(E_min, E_max, E_num)

# sweep over energies
T_num, R_num = sweep_energies(E_list, x, dx, V, dt_dimless,
                              x0=-14.0, sigma=3.0,
                              absorb_frac=absorb_frac, absorb_pow=absorb_pow, t_safety=1.6)

# analytic result
T_ana = T_analytic_rect_barrier(E_list, V0_eV=V0_eV, a_m=a_m)

# plot comparison
plt.figure(figsize=(7,4))
plt.plot(E_list, T_num, 'o', ms=3, label='numerical T_num(E)')
plt.plot(E_list, T_ana, '-', label='analytic T_ana(E)')
plt.xlabel('Energy (eV)')
plt.ylabel('Transmission')
plt.ylim(-0.02, 1.02)
plt.legend()
plt.tight_layout()
if save_fig: plt.savefig('T_vs_E_compare.png', dpi=180)
plt.show()

# optional: check if T+R is close to 1 (will be slightly less with absorbing layers)
plt.figure(figsize=(6,3.6))
plt.plot(E_list, T_num+R_num, '.', label='T_num + R_num')
plt.axhline(1.0, linestyle='--', label='1.0')
plt.xlabel('Energy (eV)')
plt.ylabel('T+R')
plt.legend()
plt.tight_layout()
if save_fig: plt.savefig('T_plus_R.png', dpi=160)
plt.show()
