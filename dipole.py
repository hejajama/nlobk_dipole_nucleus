"""
Dipole-nucleus amplitude N(x, r, b) interpolated from BK evolution output files.

Each file <directory>/Au_b_<b>.dat (b in GeV^-1) is in the nlobk solver output format:
    # comment lines
    ###minr
    ###rmultiplier        (r_i = minr * rmultiplier^i)
    ###number of r points
    ###x0
    ###y_0
    N(r_0)
    ...
    N(r_{n-1})
    ###y_1
    ...
where y is the evolution rapidity, x = x0 * exp(-y).
"""

import glob
import os
import re
import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator

_NMIN = 1e-300  # floor for N before taking the logarithm


def read_bk_file(filename):
    """
    Read one nlobk output file.

    Returns (r, y, x0, N) where N has shape (len(y), len(r)).
    """
    with open(filename) as f:
        blocks = f.read().split("###")[1:]  # drop the comment header

    minr = float(blocks[0].split()[0])
    rmult = float(blocks[1].split()[0])
    nr = int(float(blocks[2].split()[0]))
    x0 = float(blocks[3].split()[0])

    y = []
    N = []
    for block in blocks[4:]:
        vals = np.array(block.split(), dtype=float)
        if len(vals) != nr + 1:
            raise ValueError(f"{filename}: expected {nr} r points in a y block, got {len(vals) - 1}")
        y.append(vals[0])
        N.append(vals[1:])

    r = minr * rmult ** np.arange(nr)
    return r, np.array(y), x0, np.array(N)


class NucleusDipole:
    """
    Interpolator for the dipole-nucleus amplitude N(x, r, b).

    Interpolation is done on a regular grid in (ln x, ln r, b). Linear in ln x
    means the interpolation in x is logarithmic. The interpolated quantity is
    ln N (if log_amplitude=True), which is accurate for the steep power-law
    small-r behaviour N ~ r^(2 gamma); interpolating N itself overshoots there.

    Args:
        directory: directory containing files Au_b_<b>.dat
        prefix: file name prefix (default "Au_b_")
        method: interpolation method for RegularGridInterpolator
                ("linear", "cubic", ...)
        log_amplitude: interpolate ln N instead of N (default True)
    """

    def __init__(self, directory, prefix="Au_b_", method="linear", log_amplitude=True):
        pattern = re.compile(re.escape(prefix) + r"([0-9.eE+-]+)\.dat$")
        files = []
        for fname in glob.glob(os.path.join(directory, prefix + "*.dat")):
            m = pattern.search(os.path.basename(fname))
            if m:
                files.append((float(m.group(1)), fname))
        if not files:
            raise FileNotFoundError(f"No {prefix}*.dat files found in {directory}")
        files.sort()

        b_values = []
        data = []
        r_ref = y_ref = x0_ref = None
        for b, fname in files:
            r, y, x0, N = read_bk_file(fname)
            if r_ref is None:
                r_ref, y_ref, x0_ref = r, y, x0
            else:
                # Runs may stop at different final y; compare on the common y range.
                ny = min(len(y), len(y_ref))
                if not (r.shape == r_ref.shape and np.allclose(r, r_ref)
                        and np.allclose(y[:ny], y_ref[:ny]) and x0 == x0_ref):
                    raise ValueError(f"{fname}: (r, y, x0) grid differs from the other files")
                y_ref = y_ref[:ny]
            b_values.append(b)
            data.append(N)

        # Truncate all files to the common rapidity range
        ny = len(y_ref)
        if any(len(N) != ny for N in data):
            warnings.warn(f"{directory}: files have different numbers of y steps; "
                          f"truncating all to the common range y <= {y_ref[-1]:.4g}")
        data = [N[:ny] for N in data]

        self.b = np.array(b_values)
        self.r = r_ref
        self.x0 = x0_ref
        self.y = y_ref
        self.x = x0_ref * np.exp(-y_ref)
        # Grid axes must be increasing: ln x = ln x0 - y decreases with y, so flip.
        lnx = np.log(self.x)[::-1]
        N = np.array(data)  # shape (nb, ny, nr)
        N = np.transpose(N, (1, 2, 0))[::-1]  # shape (nx, nr, nb), x increasing
        self.log_amplitude = log_amplitude
        if log_amplitude:
            N = np.log(np.maximum(N, _NMIN))

        self.x_range = (self.x.min(), self.x.max())
        self.r_range = (self.r.min(), self.r.max())
        self.b_range = (self.b.min(), self.b.max())

        self._interp = RegularGridInterpolator(
            (lnx, np.log(self.r), self.b), N, method=method,
            bounds_error=False, fill_value=np.nan,
        )

    def N(self, x, r, b):
        """
        Dipole-nucleus amplitude.

        Args:
            x: Bjorken x (scalar or array)
            r: dipole size in GeV^-1 (scalar or array)
            b: impact parameter in GeV^-1 (scalar or array)

        Arrays are broadcast against each other. Returns NaN outside the
        tabulated range in x, r or b.
        """
        x, r, b = np.broadcast_arrays(np.asarray(x, dtype=float),
                                      np.asarray(r, dtype=float),
                                      np.asarray(b, dtype=float))
        pts = np.stack([np.log(x).ravel(), np.log(r).ravel(), b.ravel()], axis=-1)
        res = self._interp(pts).reshape(x.shape)
        if self.log_amplitude:
            res = np.exp(res)
        return res.item() if res.ndim == 0 else res

    __call__ = N


class ProtonDipole:
    """
    Interpolator for the (b-independent) dipole-proton amplitude N(x, r).

    Interpolation is done on a regular grid in (ln x, ln r), i.e. logarithmic in x.
    As for NucleusDipole, ln N is interpolated if log_amplitude=True.

    Args:
        directory: directory containing the proton file
        filename: file name (default "proton.dat")
        method: interpolation method for RegularGridInterpolator
        log_amplitude: interpolate ln N instead of N (default True)
    """

    def __init__(self, directory, filename="proton.dat", method="linear", log_amplitude=True):
        r, y, x0, N = read_bk_file(os.path.join(directory, filename))
        self.r = r
        self.y = y
        self.x0 = x0
        self.x = x0 * np.exp(-y)

        self.x_range = (self.x.min(), self.x.max())
        self.r_range = (self.r.min(), self.r.max())

        self.log_amplitude = log_amplitude
        if log_amplitude:
            N = np.log(np.maximum(N, _NMIN))

        # Flip so that ln x is increasing
        self._interp = RegularGridInterpolator(
            (np.log(self.x)[::-1], np.log(self.r)), N[::-1], method=method,
            bounds_error=False, fill_value=np.nan,
        )

    def N(self, x, r):
        """
        Dipole-proton amplitude.

        Args:
            x: Bjorken x (scalar or array)
            r: dipole size in GeV^-1 (scalar or array)

        Arrays are broadcast against each other. Returns NaN outside the
        tabulated range in x or r.
        """
        x, r = np.broadcast_arrays(np.asarray(x, dtype=float),
                                   np.asarray(r, dtype=float))
        pts = np.stack([np.log(x).ravel(), np.log(r).ravel()], axis=-1)
        res = self._interp(pts).reshape(x.shape)
        if self.log_amplitude:
            res = np.exp(res)
        return res.item() if res.ndim == 0 else res

    __call__ = N


if __name__ == "__main__":
    dip = NucleusDipole(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kcbk_bal", "Au"))
    print("x range", dip.x_range, "r range", dip.r_range, "b range", dip.b_range)
    for x in [1e-2, 1e-3, 1e-4]:
        print(x, [round(dip.N(x, r, 0.0), 5) for r in [0.1, 1.0, 5.0]])
