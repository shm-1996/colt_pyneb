"""
emission.py — LineSpec definition and GridEmission class.

GridEmission computes pyneb-based emission line luminosities on a simulation
grid by interpolating emissivity tables pre-built on a (Te, ne) mesh.  This
avoids calling pyneb once per cell and scales to grids with O(10^8) cells.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass(frozen=True)
class LineSpec:
    """
    Specification for a single emission line.

    Attributes
    ----------
    atom    : pyneb Atom or RecAtom instance
    lev_i   : upper level index (pyneb convention)
    lev_j   : lower level index
    ion_key : key in the ion_densities dict passed to GridEmission, e.g. 'OIII'
    wav     : reference wavelength [Å] — used for display only
    rec     : True for recombination lines (pyneb RecAtom); False for CELs
    """
    atom:    object
    lev_i:   int | None   # None for non-hydrogenic RecAtom lines (wav is used instead)
    lev_j:   int | None
    ion_key: str
    wav:     float        # reference wavelength [Å]; also passed to getEmissivity for non-hydrogenic rec lines
    rec:     bool = False


class GridEmission:
    """
    Compute emission line luminosities on a simulation grid using pyneb.

    Emissivity coefficients ε(Te, ne) [erg s⁻¹ cm³] are pre-tabulated on a
    log-spaced (Te, ne) grid and bi-linearly interpolated.  Luminosity per cell:

        L_cell = ε(Te, ne) × n_ion × n_e × V_cell   [erg/s]

    Parameters
    ----------
    Te : (N,) array
        Electron temperature [K] for each grid cell.
    ne : (N,) array
        Electron number density [cm⁻³] for each grid cell.
    ion_densities : dict[str, (N,) array]
        Ion number densities [cm⁻³], keyed by species label (e.g. 'OIII', 'HII').
        Must include an entry for every `ion_key` used by the requested lines.
    vol_cell : float or (N,) array
        Cell volume [cm³].  A scalar assumes uniform cells.
    Te_range : (T_min, T_max)
        Bounds of the temperature interpolation grid [K].  Default (3e3, 1e5).
    ne_range : (n_min, n_max)
        Bounds of the density interpolation grid [cm⁻³].  Default (1e-2, 1e7).
    n_Te, n_ne : int
        Number of grid points along each axis (logarithmically spaced).
    """

    def __init__(self, Te, ne, ion_densities: dict, vol_cell,
                 Te_range: tuple = (3e3, 1e5),
                 ne_range: tuple = (1e-2, 1e7),
                 n_Te: int = 80, n_ne: int = 80):
        self.Te   = np.asarray(Te,  dtype=float)
        self.ne   = np.asarray(ne,  dtype=float)
        self.ions = {k: np.asarray(v, dtype=float) for k, v in ion_densities.items()}
        self.vol  = vol_cell
        self._log_Te = np.linspace(np.log10(Te_range[0]), np.log10(Te_range[1]), n_Te)
        self._log_ne = np.linspace(np.log10(ne_range[0]), np.log10(ne_range[1]), n_ne)
        self._cache: dict[LineSpec, RegularGridInterpolator] = {}

    # ------------------------------------------------------------------
    # Emissivity table
    # ------------------------------------------------------------------

    def _get_interp(self, spec: LineSpec) -> RegularGridInterpolator:
        """Return a cached log-log 2D interpolator for spec's emissivity."""
        if spec in self._cache:
            return self._cache[spec]

        Te1d = 10.0 ** self._log_Te
        ne1d = 10.0 ** self._log_ne
        # pyneb.getEmissivity(tem_array, den_array, ...) returns the outer product
        # as a 2D array of shape (len(tem), len(den)) — no meshgrid needed.
        # Non-hydrogenic RecAtom lines (lev_i is None) are identified by wavelength.
        if spec.rec and spec.lev_i is None:
            eps = np.array(
                spec.atom.getEmissivity(Te1d, ne1d, wave=spec.wav),
                dtype=float,
            )
        else:
            eps = np.array(
                spec.atom.getEmissivity(Te1d, ne1d, spec.lev_i, spec.lev_j),
                dtype=float,
            )
        log_eps = np.log10(np.abs(eps) + 1e-300)

        # Some pyneb tables (e.g. SH95 H recombination, valid only to ~30 000 K)
        # return NaN outside their Te range.  Fill with the nearest valid value
        # so the interpolator extrapolates from the boundary rather than
        # propagating NaN into cells that lie outside the table coverage.
        if np.any(np.isnan(log_eps)):
            from scipy.ndimage import distance_transform_edt
            idx = distance_transform_edt(
                np.isnan(log_eps), return_distances=False, return_indices=True
            )
            log_eps = log_eps[tuple(idx)]

        interp = RegularGridInterpolator(
            (self._log_Te, self._log_ne),
            log_eps,
            method='linear', bounds_error=False, fill_value=None,
        )
        self._cache[spec] = interp
        return interp

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def cell_luminosities(self, lines, mask=None) -> dict[str, np.ndarray]:
        """
        Compute per-cell emission line luminosities [erg/s].

        Parameters
        ----------
        lines : list[str] | dict[str, LineSpec]
            Either a list of LINE_LIBRARY keys (e.g. ['OIII_5007', 'Ha'])
            or a custom {name: LineSpec} dict for lines not in the library.
        mask : array-like of bool or int, optional
            Restrict computation to a subset of cells (e.g. ionized region).
            Cells outside the mask are set to zero in the output arrays.

        Returns
        -------
        dict[str, (N,) ndarray]  — per-cell luminosity [erg/s]
        """
        line_specs = self._resolve_lines(lines)
        N = len(self.Te)

        if mask is None:
            idx = np.arange(N)
        else:
            mask = np.asarray(mask)
            idx  = np.where(mask)[0] if mask.dtype == bool else mask.astype(int)

        # Query points in log space, guarded against zero/negative inputs
        log_Te = np.log10(np.maximum(self.Te[idx], 1.0))
        log_ne = np.log10(np.maximum(self.ne[idx], 1e-300))
        pts    = np.column_stack([log_Te, log_ne])

        print("Max temperature and density = ",np.max(self.Te[idx]),np.max(self.ne[idx]))

        chunk  = 5_000_000
        result = {}

        for name, spec in line_specs.items():
            if spec.ion_key not in self.ions:
                raise KeyError(
                    f"ion_key '{spec.ion_key}' (required by line '{name}') not found in "
                    f"ion_densities.  Available keys: {sorted(self.ions)}"
                )
            interp = self._get_interp(spec)
            L      = np.zeros(N, dtype=float)
            n_ion  = self.ions[spec.ion_key]
            scalar_vol = np.ndim(self.vol) == 0

            for i in range(0, len(idx), chunk):
                sl  = slice(i, min(i + chunk, len(idx)))
                ci  = idx[sl]
                eps = 10.0 ** interp(pts[sl])
                vol = self.vol if scalar_vol else self.vol[ci]
                L[ci] = eps * n_ion[ci] * self.ne[ci] * vol

            result[name] = L

        return result

    def integrated(self, lines, mask=None) -> dict[str, float]:
        """
        Sum cell luminosities to give total line luminosities [erg/s].

        Parameters
        ----------
        lines : list[str] | dict[str, LineSpec]
            Lines to compute (same as cell_luminosities).
        mask : array-like, optional
            Same masking as cell_luminosities; only masked cells contribute.

        Returns
        -------
        dict[str, float]  — total luminosity per line [erg/s]
        """
        return {name: float(arr.sum())
                for name, arr in self.cell_luminosities(lines, mask=mask).items()}

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _resolve_lines(self, lines) -> dict[str, LineSpec]:
        """Convert a list of library keys or a custom dict to {name: LineSpec}."""
        if isinstance(lines, dict):
            return lines
        from .line_library import LINE_LIBRARY
        missing = [k for k in lines if k not in LINE_LIBRARY]
        if missing:
            raise KeyError(f"Lines not found in LINE_LIBRARY: {missing}.  "
                           f"Pass a {{name: LineSpec}} dict to use custom lines.")
        return {k: LINE_LIBRARY[k] for k in lines}
