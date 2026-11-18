#!/usr/bin/env python3
"""
Dipole-nucleus initial conditions for (KC)BK evolution at fixed impact parameter,
Eq. (19) of arXiv:2311.10491 (optical Glauber):

  N_A(r,b) = 1 - exp[ -A T_A(b) sigma0/2 * (r^2 Qs0^2)^gamma / 4 * ln(1/(r Lambda) + e_c e) ]

used at all b (the large-b prescription N_A = A T_A sigma0/2 * N_p of the paper is NOT applied).

Parameters: KCBK + PD fit of arXiv:2506.00487, MAP values (Table 1):
  Qs0^2 = 0.090 GeV^2, gamma = 1.6, e_c = 1, C^2 = 663, sigma0/2 = 20.7 mb, x0 = 1
(C^2 enters only the running coupling, i.e. the nlobk command line.)

Output format readable by nlobk `-ic FILE`: lines "r[GeV^-1] N".
"""
import numpy as np
from scipy.integrate import quad

A = 197                      # gold
QS0SQR = 0.090               # GeV^2
GAMMA = 1.6
EC = 1.0
SIGMA02_MB = 20.7            # mb
LAMBDA_QCD = 0.241           # GeV, same as nlobk config::LAMBDAQCD
FMGEV = 0.1973269804         # GeV fm
MB_TO_GEVM2 = 0.1 / FMGEV**2 # 1 mb = 0.1 fm^2

SIGMA02 = SIGMA02_MB * MB_TO_GEVM2   # GeV^-2

# Woods-Saxon (as in arXiv:2311.10491)
D_WS = 0.54 / FMGEV                                        # GeV^-1
R_WS = (1.12 * A**(1/3) - 0.86 * A**(-1/3)) / FMGEV        # GeV^-1


def T_A_unnorm(b):
    f = lambda z: 1.0 / (1.0 + np.exp((np.sqrt(b*b + z*z) - R_WS) / D_WS))
    return 2.0 * quad(f, 0, R_WS + 30*D_WS + b, limit=200)[0]


NORM = 2*np.pi * quad(lambda b: b * T_A_unnorm(b), 0, R_WS + 30*D_WS, limit=200)[0]


def T_A(b):
    """Thickness function in GeV^2 (inverse area), normalised to int d^2b T_A = 1 (b in GeV^-1)."""
    return T_A_unnorm(b) / NORM


def exponent_p(r):
    """Proton exponent of Eq. (10)."""
    return (r*r*QS0SQR)**GAMMA / 4.0 * np.log(1.0/(r*LAMBDA_QCD) + EC*np.e)


def N_A(r, b):
    return -np.expm1(-A * T_A(b) * SIGMA02 * exponent_p(r))


rvals = np.logspace(-7, 3, 600)   # GeV^-1; nlobk grid is 1e-5 ... 30 by default

for b in range(0, 33, 2):
    c = A * T_A(b) * SIGMA02
    nvals = N_A(rvals, b)
    mode = "Eq.(19), optical Glauber (used at all b)"
    fname = f"./kcbk_pd/Au/ic_Au_b_{b}.dat"
    with open(fname, "w") as f:
        f.write(f"# Dipole-nucleus IC, Eq. (19) of arXiv:2311.10491, A={A} (Au), b={b} GeV^-1 = {b*FMGEV:.4f} fm\n")
        f.write(f"# KCBK+PD MAP (arXiv:2506.00487): Qs0^2={QS0SQR} GeV^2 gamma={GAMMA} e_c={EC} "
                f"sigma0/2={SIGMA02_MB} mb={SIGMA02:.4f} GeV^-2 C^2=663 x0=1.00 Lambda_QCD={LAMBDA_QCD} GeV\n")
        f.write(f"# Woods-Saxon R_A={R_WS*FMGEV:.4f} fm d=0.54 fm; T_A(b)={T_A(b):.6e} GeV^2 = {T_A(b)/FMGEV**2:.6e} fm^-2; "
                f"A T_A(b) sigma0/2 = {c:.6f}\n")
        f.write(f"# Mode: {mode}\n")
        f.write("# r [GeV^-1]  N(r)\n")
        for r, n in zip(rvals, nvals):
            f.write(f"{r:.10e} {n:.10e}\n")
    print(f"b={b:2d} GeV^-1  T_A={T_A(b):.4e} GeV^2  A T_A sigma0/2={c:8.4f}  {mode.split(':')[0]}")
