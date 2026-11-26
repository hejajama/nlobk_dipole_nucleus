# BK-evolved dipole amplitudes for protons and nuclei

This repository contains dipole–proton and dipole–nucleus scattering amplitudes
obtained by solving the Balitsky–Kovchegov (BK) equation.

- **Proton:** the dipole–proton amplitude $N_p(r, x)$ does not depend on impact parameter.
- **Nucleus:** the dipole–nucleus amplitude $N_A(r, x, b)$ is computed by evolving the BK
  equation separately at each impact parameter $b$. Different impact parameters do not
  mix during the evolution.

The evolution is computed with the [nlobk](https://github.com/hejajama/nlobk) solver.

## Available dipoles

| Directory   | Description | Contents |
|-------------|-------------|----------|
| `kcbk_bal/` | KCBK, Balitsky + smallest-dipole running coupling, MAP parameters | `proton.dat`; `Au/`, `Pb/`: `<A>_b_<b>.dat` (evolved), `ic_<A>_b_<b>.dat` (initial conditions) |
| `kcbk_pd/`  | KCBK, parent-dipole running coupling, MAP parameters | `proton.dat`; `Au/`, `Pb/`: `<A>_b_<b>.dat` (evolved), `ic_<A>_b_<b>.dat` (initial conditions) |
| `nlobk_mv_balsd/` | NLO BK, MV initial condition, Balitsky + smallest-dipole running coupling | `proton.dat`; `Au/`, `Pb/`: `<A>_b_<b>.dat` (evolved), `ic_<A>_b_<b>.dat` (initial conditions) |

Each fit directory holds the proton dipole `proton.dat` and one subdirectory per nucleus,
named after the nucleus (`Au/`, `Pb/`).

The KCBK dipoles are taken from the fit to HERA data in
C. Casuga, H. Hänninen and H. Mäntysaari, Phys. Rev. D 112, 034003 (2025),
[arXiv:2506.00487](https://arxiv.org/abs/2506.00487), using the MAP (maximum a
posteriori) parameter values. KCBK means the BK equation with a kinematical constraint
from [arXiv:1401.0313](https://arxiv.org/abs/1401.0313).

The NLO BK dipoles (`nlobk_mv_balsd`) use the fit to HERA data in
[arXiv:2604.22332](https://arxiv.org/abs/2604.22332): MV initial condition
($\gamma = 1$), NLO BK with resummation of double and single transverse logarithms
($K_\mathrm{sub} = 0.65$), Balitsky running coupling in the LO kernel and
smallest-dipole running coupling in the NLO terms, $N_f = 4$ in the $n_f K_f$ term of the NLO BK,
$r$-dependent $n_f$ when evaluating $\alpha_s$.

Dipoles from other setups will be added later.

### Grids

- Dipole size: logarithmic grid, same for the nucleus and the proton:
  $r \in [10^{-5}, 30]\ \mathrm{GeV}^{-1}$ (`kcbk_bal`, 101 points),
  $r \in [10^{-5}, 200]\ \mathrm{GeV}^{-1}$ (`kcbk_pd`, 101 points) and
  $r \in [2\cdot 10^{-5}, 150]\ \mathrm{GeV}^{-1}$ (`nlobk_mv_balsd`, 121 points).
- Impact parameter (nucleus): $b = 0, 2, 4, \ldots, 32\ \mathrm{GeV}^{-1}$, with one file per $b$.
- Evolution rapidity: $y = \ln(x_0/x)$ from 0 to 15.

All dimensionful quantities are in GeV units, so $r$ and $b$ are given in $\mathrm{GeV}^{-1}$.

## Initial conditions

The proton initial condition at $x = x_0$ is the MV$^\gamma$ form

$$
N_p(r) = 1 - \exp\left[ -\frac{(r^2 Q_{s,0}^2)^\gamma}{4}
\ln\left( \frac{1}{r\Lambda_\mathrm{QCD}} + e_c \cdot e \right) \right].
$$

The dipole–nucleus initial condition uses the optical Glauber model, following
T. Lappi and H. Mäntysaari, Phys. Rev. D 88, 114020 (2013),
[arXiv:1309.6963](https://arxiv.org/abs/1309.6963):

$$
N_A(r, b) = 1 - \exp\left[ -A\, T_A(b)\, \frac{\sigma_0}{2}\,
\frac{(r^2 Q_{s,0}^2)^\gamma}{4}
\ln\left( \frac{1}{r\Lambda_\mathrm{QCD}} + e_c \cdot e \right) \right].
$$

Here $\sigma_0/2$ is the proton transverse area from the fit. $T_A(b)$ is the nuclear
thickness function, obtained from a Woods–Saxon density and normalized as
$\int \mathrm{d}^2 b\, T_A(b) = 1$. This form is used at every impact parameter, and the
script `generate_ic.py` generates the initial-condition files. Choose the nucleus and the fit
on the command line:

```bash
python3 generate_ic.py --nucleus Pb --fit pd    # writes kcbk_pd/Pb/ic_Pb_b_<b>.dat
python3 generate_ic.py --nucleus Au --fit bal   # writes kcbk_bal/Au/ic_Au_b_<b>.dat
```

### Parameters

The values for $Q_{s,0}^2$, $\gamma$, $e_c$, $\sigma_0/2$ and $C^2$ are the MAP values from
[arXiv:2506.00487](https://arxiv.org/abs/2506.00487).

For `nlobk_mv_balsd` they are taken from
[arXiv:2604.22332](https://arxiv.org/abs/2604.22332).

| Parameter | `kcbk_bal` | `kcbk_pd` | `nlobk_mv_balsd` |
|-----------|-----------:|----------:|-----------------:|
| $Q_{s,0}^2$ [GeV²] | 0.124 | 0.090 | 0.00721 |
| $\gamma$ | 3.23 | 1.6 | 1 |
| $e_c$ | 1 | 1 | 1 |
| $\sigma_0/2$ [mb] | 9.08 | 20.7 | 22.8 |
| $C^2$ (running coupling) | 1.74 | 663 | TODO |
| $x_0$ | 1 | 1 | 1 |
| $\Lambda_\mathrm{QCD}$ [GeV] | 0.241 | 0.241 | 0.241 |
| $N_f$ | 3 | 3 | 4 |
| BK equation | LO, KCBK | LO, KCBK | NLO, double + single log resummed |
| Running-coupling prescription | Balitsky (LO) | parent dipole | Balitsky (LO kernel), smallest dipole (NLO kernel) |

Nuclei:

| Parameter | Au | Pb |
|-----------|---:|---:|
| $A$ | 197 | 208 |
| Woods–Saxon radius $R_A = (1.12 A^{1/3} - 0.86 A^{-1/3})$ fm | 6.37 fm | 6.49 fm |
| Woods–Saxon surface thickness $d$ | 0.54 fm | 0.54 fm |

### Evolution

Each impact parameter is evolved separately with [nlobk](https://github.com/hejajama/nlobk)
up to $y = 15$, using the LO kernel with the kinematical constraint:

```bash
nlobk -ic FILE kcbk_bal/Pb/ic_Pb_b_<b>.dat -order lo -rc balitsky -alphas_scaling 1.74 -nf 3 -maxy 15 -output kcbk_bal/Pb/Pb_b_<b>.dat
nlobk -ic FILE kcbk_pd/Pb/ic_Pb_b_<b>.dat  -order lo -rc parent   -alphas_scaling 663  -nf 3 -maxy 15 -output kcbk_pd/Pb/Pb_b_<b>.dat
```

The $r$ grid is fixed when nlobk is compiled (`MAXR` in `src/nlobk_config.cpp`):
`kcbk_bal` uses the default `MAXR = 30` and `kcbk_pd` uses `MAXR = 200`.
nlobk writes `###0.01` as the $x_0$ line; it is set to `###1.00` in the files here.

### Running coupling

The coupling runs in coordinate space:

$$
\alpha_s(r^2) = \frac{4\pi}{\beta_0 \ln\left\{ \left[ \left(\frac{\mu_0^2}{\Lambda_\mathrm{QCD}^2}\right)^{1/c}
+ \left( \frac{4 C^2}{r^2 \Lambda_\mathrm{QCD}^2} \right)^{1/c} \right]^{c} \right\}},
\qquad \beta_0 = \frac{11 N_c - 2 N_f}{3}.
$$

Here $\mu_0/\Lambda_\mathrm{QCD} = 2.5$ and $c = 0.2$. These two values make the coupling
freeze smoothly in the infrared at $\alpha_s \approx 0.76$ for $N_f = 3$. The fit parameter
$C^2$ sets the relation between the transverse coordinate $r$ and the momentum scale.

## File format

Each `.dat` file uses the nlobk output format. It starts with `#` comment lines that
describe the setup, followed by blocks separated by `###`:

```
###minr            r_i = minr * rmultiplier^i
###rmultiplier
###number of r points
###x0
###y_0
N(r_0)
...
N(r_{n-1})
###y_1
...
```

## Example

`dipole.py` provides interpolators for the amplitudes: `ProtonDipole` and
`NucleusDipole`. Interpolation is linear in $(\ln x, \ln r, b)$ and acts on $\ln N$.
Outside the tabulated range, the interpolators return `NaN`.

```python
import numpy as np
from dipole import ProtonDipole, NucleusDipole

proton = ProtonDipole("kcbk_bal")       # reads kcbk_bal/proton.dat
nucleus = NucleusDipole("kcbk_bal", nucleus="Au")  # reads kcbk_bal/Au/Au_b_*.dat
# lead = NucleusDipole("kcbk_bal", nucleus="Pb")  # reads kcbk_bal/Pb/Pb_b_*.dat

# The tables are stored on a grid of evolution rapidity y, measured from the
# initial condition. dip.x0 is the reference value written in the data file,
# so x = dip.x0 * exp(-y).
y = 4.0
r = np.array([0.1, 1.0, 5.0])     # dipole size [GeV^-1]

# Dipole-proton amplitude N_p(r) at rapidity y
Np = proton.N(proton.x0 * np.exp(-y), r)

# Dipole-nucleus amplitude N_A(r, b) at rapidity y and b = 0, 10 GeV^-1
NA_central = nucleus.N(nucleus.x0 * np.exp(-y), r, 0.0)
NA_b10 = nucleus.N(nucleus.x0 * np.exp(-y), r, 10.0)

print("r [GeV^-1]  :", r)
print("N_p         :", Np)
print("N_A(b=0)    :", NA_central)
print("N_A(b=10)   :", NA_b10)

# The inputs broadcast against each other, so N_A(r, b) on a 2D grid is:
R, B = np.meshgrid(np.logspace(-2, 1, 50), np.linspace(0, 32, 33), indexing="ij")
NA_grid = nucleus.N(nucleus.x0 * np.exp(-y), R, B)  # shape (50, 33)
```

If the evolved files for different $b$ reach different final rapidities,
`NucleusDipole` truncates all of them to the common range and issues a warning.

`verify.ipynb` contains further checks, for example a comparison of
$\int \mathrm{d}^2 b\, N_A$ with $A\,\frac{\sigma_0}{2} N_p$ in the dilute limit.

## References

- KCBK fit and MAP parameters: C. Casuga, H. Hänninen, H. Mäntysaari,
  [arXiv:2506.00487](https://arxiv.org/abs/2506.00487)
- NLO BK fit with MV initial condition (`nlobk_mv_balsd`):
  [arXiv:2604.22332](https://arxiv.org/abs/2604.22332)
- Optical Glauber initial condition for nuclei: T. Lappi, H. Mäntysaari,
  [arXiv:1309.6963](https://arxiv.org/abs/1309.6963)
- BK solver: [github.com/hejajama/nlobk](https://github.com/hejajama/nlobk)
