"""
colt.py — COLT/FLASH snapshot loader for use with nebular_em.

COLTSnapshot reads the states-teq_NNNN.hdf5 and colt_NNNN.hdf5 files produced
by a COLT radiative transfer run and exposes the gas grid arrays and ion number
densities that GridEmission and Diagnostics need.

Solar abundances follow Grevesse et al. (2010) and match the normalisation
used in COLT.  Ion fractions in the states file are stored as fractions of the
parent-element number density (y_X = n_X / n_element_total) and are converted
to absolute number densities [cm⁻³] here.

Other COLT-specific functionality (e.g. loading COLT-output emission line
luminosities to feed directly into Diagnostics) will be added to this module.
"""

import os
import numpy as np
import h5py

PROTONMASS: float = 1.6726e-24   # g
_PC:        float = 3.085677581e18  # cm per parsec

# Solar mass fractions Z_sun,X — Grevesse et al. (2010), same values as COLT
SOLAR_MASS_FRAC: dict[str, float] = {
    'C':  0.00236729,
    'N':  0.000693438,
    'O':  0.00573805,
    'Ne': 0.00125773,
    'Mg': 0.000708545,
    'Si': 0.000665509,
    'S':  0.00030953,
    'Fe': 0.00129317,
}

# Atomic masses [amu]
ATOMIC_MASS: dict[str, float] = {
    'H':  1.008,  'He': 4.003,  'C':  12.011, 'N':  14.007,
    'O':  15.999, 'Ne': 20.180, 'Mg': 24.305, 'Si': 28.085,
    'S':  32.060, 'Fe': 55.845,
}


def solar_number_abundance(element: str, X_H: float = 0.76) -> float:
    """
    Return n_X / n_H at solar metallicity using Grevesse et al. (2010).

    From the mass-fraction definition  Z_sun,X = X_H * (m_X/m_H) * (n_X/n_H):
        n_X/n_H = Z_sun,X * m_H / (X_H * m_X)
    """
    return SOLAR_MASS_FRAC[element] * ATOMIC_MASS['H'] / (X_H * ATOMIC_MASS[element])


# Maps nebular_em ion_key -> (HDF5 dataset name in states file, element symbol)
#
# Abundance convention:
#   'H'  : n_ion = nH * x                             (x is fraction of nH)
#   'He' : n_ion = nH * x * X_He / (4 * X_H)         (x is fraction of total He)
#   else : n_ion = nH * x * (n_X/nH)_solar * Zgas_Zsun
ION_MAP: dict[str, tuple[str, str]] = {
    # Hydrogen
    'HI'   : ('x_HI',    'H'),
    'HII'  : ('x_HII',   'H'),
    # Helium
    'HeI'  : ('x_HeI',   'He'),
    'HeII' : ('x_HeII',  'He'),
    # Carbon
    'CI'   : ('x_CI',    'C'),
    'CII'  : ('x_CII',   'C'),
    'CIII' : ('x_CIII',  'C'),
    'CIV'  : ('x_CIV',   'C'),
    # Nitrogen
    'NI'   : ('x_NI',    'N'),
    'NII'  : ('x_NII',   'N'),
    'NIII' : ('x_NIII',  'N'),
    'NIV'  : ('x_NIV',   'N'),
    'NV'   : ('x_NV',    'N'),
    # Oxygen
    'OI'   : ('x_OI',    'O'),
    'OII'  : ('x_OII',   'O'),
    'OIII' : ('x_OIII',  'O'),
    'OIV'  : ('x_OIV',   'O'),
    # Neon
    'NeI'  : ('x_NeI',   'Ne'),
    'NeII' : ('x_NeII',  'Ne'),
    'NeIII': ('x_NeIII', 'Ne'),
    'NeIV' : ('x_NeIV',  'Ne'),
    # Magnesium
    'MgI'  : ('x_MgI',   'Mg'),
    'MgII' : ('x_MgII',  'Mg'),
    'MgIII': ('x_MgIII', 'Mg'),
    # Silicon
    'SiI'  : ('x_SiI',   'Si'),
    'SiII' : ('x_SiII',  'Si'),
    'SiIII': ('x_SiIII', 'Si'),
    'SiIV' : ('x_SiIV',  'Si'),
    # Sulfur
    'SI'   : ('x_SI',    'S'),
    'SII'  : ('x_SII',   'S'),
    'SIII' : ('x_SIII',  'S'),
    'SIV'  : ('x_SIV',   'S'),
    'SV'   : ('x_SV',    'S'),
    # Iron
    'FeI'  : ('x_FeI',   'Fe'),
    'FeII' : ('x_FeII',  'Fe'),
    'FeIII': ('x_FeIII', 'Fe'),
    'FeIV' : ('x_FeIV',  'Fe'),
    'FeV'  : ('x_FeV',   'Fe'),
    'FeVI' : ('x_FeVI',  'Fe'),
}


class COLTSnapshot:
    """
    Load a single COLT/FLASH snapshot and expose grid arrays and ion densities.

    Reads states-teq_NNNN.hdf5 (ionization equilibrium state) and
    colt_NNNN.hdf5 (density, velocity, stellar particles).  Ion number
    densities are computed on demand and cached.

    Parameters
    ----------
    snap_dir   : directory containing the snapshot HDF5 files
    snapnum    : snapshot number (integer), used to build filenames
    Zgas_Zsun  : uniform gas metallicity in solar units — must match the COLT run
    N          : grid cells per side (default 512 for 512³ runs)
    X_H        : hydrogen mass fraction assumed in the simulation (default 0.76)
    thermal_eq : if True (default) read states-teq_NNNN.hdf5; else states_NNNN.hdf5

    Examples
    --------
    >>> snap_dir = '/path/to/colt_dir'
    >>> with COLTSnapshot(snap_dir, snapnum=150, Zgas_Zsun=0.1) as snap:
    ...     grid = snap.to_grid_emission(['OIII', 'OII', 'HII', 'SII'])
    ...     L    = grid.integrated(['OIII_5007', 'OIII_4363', 'OII_3726', 'OII_3729'],
    ...                             mask=snap.ionized_mask())
    ...     from nebular_em import Diagnostics
    ...     Te, ne = Diagnostics(L).iterate_Te_ne(
    ...                 ('OIII_4363', 'OIII_5007'), ('OII_3726', 'OII_3729'))
    """

    def __init__(self, snap_dir: str, snapnum: int, Zgas_Zsun: float,
                 N: int = 512, X_H: float = 0.76, thermal_eq: bool = True):
        self.snap_dir   = snap_dir
        self.snapnum    = snapnum
        self.Zgas_Zsun  = Zgas_Zsun
        self.N          = N
        self.X_H        = X_H
        self.thermal_eq = thermal_eq

        prefix      = 'states-teq' if thermal_eq else 'states'
        states_path = os.path.join(snap_dir, f'{prefix}_0{snapnum:03d}.hdf5')
        colt_path   = os.path.join(snap_dir, f'colt_0{snapnum:03d}.hdf5')

        for p in (states_path, colt_path):
            if not os.path.exists(p):
                raise FileNotFoundError(f'Snapshot file not found: {p}')

        self._states = h5py.File(states_path, 'r')
        self._colt   = h5py.File(colt_path,   'r')
        self._cache: dict = {}

    # ── Gas grid arrays ────────────────────────────────────────────────────────

    @property
    def Te(self) -> np.ndarray:
        """Equilibrium electron temperature [K], shape (N³,)."""
        return self._load('Te', lambda: self._states['T_eq'][:])

    @property
    def nH(self) -> np.ndarray:
        """Hydrogen nucleon density [cm⁻³] = nHI + nHII + 2 nH₂."""
        return self._load('nH', lambda: self._colt['rho'][:] * self.X_H / PROTONMASS)

    @property
    def ne(self) -> np.ndarray:
        """Electron number density [cm⁻³] from the post-COLT ionization state."""
        return self._load('ne', lambda: self._states['x_e'][:] * self.nH)

    @property
    def x_HII(self) -> np.ndarray:
        """HII ionization fraction = n_HII / n_H_nuclei."""
        return self._load('x_HII', lambda: self._states['x_HII'][:])

    @property
    def bbox(self) -> np.ndarray:
        """Simulation bounding box [cm], shape (2, 3): [lower, upper] corners."""
        return self._load('bbox', lambda: self._colt['bbox'][:])

    @property
    def vol_cell(self) -> float:
        """Volume of a single (uniform cubic) grid cell [cm³]."""
        return self._load('vol_cell', lambda:
                          ((self.bbox[1, 0] - self.bbox[0, 0]) / self.N) ** 3)

    # ── Ion number densities ───────────────────────────────────────────────────

    def ion_densities(self, species: list[str]) -> dict[str, np.ndarray]:
        """
        Compute number densities [cm⁻³] for the requested ion species.

        Uses Grevesse et al. (2010) solar abundances and scales all metal
        species by Zgas_Zsun.  Results are cached after the first call.

        Parameters
        ----------
        species : list of keys from ION_MAP, e.g. ['OIII', 'OII', 'HII', 'SII']

        Returns
        -------
        dict[str, (N³,) ndarray]  — number density [cm⁻³] per species
        """
        unknown = [s for s in species if s not in ION_MAP]
        if unknown:
            raise KeyError(f"Unknown species: {unknown}.  "
                           f"Supported keys: {sorted(ION_MAP)}")
        nH = self.nH
        return {key: self._load(f'n_{key}', lambda k=key: self._n_ion(k, nH))
                for key in species}

    def _n_ion(self, ion_key: str, nH: np.ndarray) -> np.ndarray:
        hdf5_key, element = ION_MAP[ion_key]
        x = self._states[hdf5_key][:]
        if element == 'H':
            return nH * x
        if element == 'He':
            # x_HeI/II are fractions of total helium; convert to fraction of nH
            return nH * x * (1.0 - self.X_H) / (4.0 * self.X_H)
        return nH * x * solar_number_abundance(element, self.X_H) * self.Zgas_Zsun

    # ── Convenience ───────────────────────────────────────────────────────────

    def ionized_mask(self, threshold: float = 0.2) -> np.ndarray:
        """Boolean mask of cells with x_HII > threshold (default 0.2)."""
        return self.x_HII > threshold

    def to_grid_emission(self, species: list[str], **kwargs):
        """
        Construct a GridEmission from this snapshot.

        Parameters
        ----------
        species  : ion_key list passed to ion_densities()
        **kwargs : forwarded to GridEmission (e.g. Te_range, ne_range)

        Returns
        -------
        nebular_em.GridEmission
        """
        from .emission import GridEmission
        return GridEmission(
            Te=self.Te,
            ne=self.ne,
            ion_densities=self.ion_densities(species),
            vol_cell=self.vol_cell,
            **kwargs,
        )

    def plot_emission_map(
        self,
        line: str,
        line2: str | None = None,
        proj_axis: int = 2,
        mask_ionized: bool = True,
        cmap: str | None = None,
        vmin=None,
        vmax=None,
        fig=None,
        ax=None,
    ):
        """
        Plot a pyneb-computed projected emission line map for this snapshot.

        Per-cell luminosities from GridEmission are summed along proj_axis to
        produce a 2-D surface map.  If line2 is given, plots log₁₀(line/line2)
        instead of the raw luminosity of line.

        Parameters
        ----------
        line : str
            LINE_LIBRARY key for the numerator (or sole) line, e.g. 'OIII_5007'.
        line2 : str, optional
            LINE_LIBRARY key for the denominator.  When provided the image shows
            log₁₀(line / line2); default colormap switches to 'RdBu_r'.
        proj_axis : int
            Axis to project along (0=x, 1=y, 2=z).  Default 2.
        mask_ionized : bool
            Restrict to ionized cells (x_HII > 0.2).  Default True.
        cmap : str, optional
            Colormap passed to imshow.  Defaults to 'magma' (single line) or
            'RdBu_r' (line ratio).
        vmin, vmax : float, optional
            Color scale limits.  For a single line these are in erg s⁻¹ and
            passed to LogNorm; for a ratio they are the log₁₀ ratio limits
            passed to Normalize (default: symmetric about zero).  When omitted
            the limits are derived automatically from the data.
        fig, ax : optional
            Existing figure and axes to draw into.  If either is None a new
            figure is created.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ax : matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from .line_library import LINE_LIBRARY, LINE_LABELS, LINE_ALIASES

        # Resolve aliases → constituent LINE_LIBRARY keys
        def _keys(name):
            return LINE_ALIASES.get(name, [name])

        def _ion_keys(name):
            return {LINE_LIBRARY[k].ion_key for k in _keys(name)}

        def _proj(name, Ls):
            return sum(Ls[k].reshape(self.N, self.N, self.N) for k in _keys(name)).sum(axis=proj_axis)

        mask   = self.ionized_mask() if mask_ionized else None
        label  = LINE_LABELS.get(line,  line)
        label2 = LINE_LABELS.get(line2, line2) if line2 is not None else None

        # Build one GridEmission covering all species needed by line (and line2)
        names   = _keys(line) + (_keys(line2) if line2 else [])
        species = list(_ion_keys(line) | (_ion_keys(line2) if line2 else set()))
        grid    = self.to_grid_emission(species)
        Ls      = grid.cell_luminosities(names, mask=mask)

        L1_2d = _proj(line, Ls)
        if line2 is None:
            L_2d = L1_2d
        else:
            L2_2d = _proj(line2, Ls)
            with np.errstate(divide='ignore', invalid='ignore'):
                L_2d = np.where((L1_2d > 0) & (L2_2d > 0),
                                np.log10(L1_2d / L2_2d), np.nan)

        ax_idx = [i for i in range(3) if i != proj_axis]
        extent = [self.bbox[0, ax_idx[0]] / _PC, self.bbox[1, ax_idx[0]] / _PC,
                  self.bbox[0, ax_idx[1]] / _PC, self.bbox[1, ax_idx[1]] / _PC]

        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))

        if line2 is None:
            _cmap = cmap or 'magma'
            pos   = L_2d[L_2d > 0]
            _vmax = vmax if vmax is not None else (float(np.percentile(pos, 95)) if len(pos) else 1.0)
            _vmin = vmin if vmin is not None else _vmax * 1e-5
            norm  = mcolors.LogNorm(vmin=_vmin, vmax=_vmax)
            cbar_label = r'erg s$^{-1}$'
        else:
            _cmap  = cmap or 'RdBu_r'
            finite = L_2d[np.isfinite(L_2d)]
            _vmax  = vmax if vmax is not None else (float(np.percentile(np.abs(finite), 95)) if len(finite) else 1.0)
            _vmin  = vmin if vmin is not None else -_vmax
            norm   = mcolors.Normalize(vmin=_vmin, vmax=_vmax)
            cbar_label = rf'$\log_{{10}}$({label} / {label2})'

        im = ax.imshow(
            L_2d.T, origin='lower', extent=extent,
            norm=norm, cmap=_cmap, interpolation='nearest',
        )
        fig.colorbar(im, ax=ax, label=cbar_label)

        ax_names = ['x', 'y', 'z']
        ax.set_xlabel(f'{ax_names[ax_idx[0]]} (pc)', fontsize=12)
        ax.set_ylabel(f'{ax_names[ax_idx[1]]} (pc)', fontsize=12)
        return fig, ax

    def plot_phase_diagram(
        self,
        weights: np.ndarray | None = None,
        bins: int | tuple[int, int] = 128,
        ne_range: tuple[float, float] | None = None,
        Te_range: tuple[float, float] | None = None,
        mask_ionized: bool = True,
        cmap: str | None = None,
        norm_log: bool = True,
        cbar_label: str | None = None,
        cbar: bool = True,
        fig=None,
        ax=None,
    ):
        """
        Plot a 2D weighted histogram in ne–Te phase space with log-scaled axes.

        Each cell contributes to the bin corresponding to its (ne, Te) with the
        given weight.  Bins are log-spaced to match the log-scaled axes.  The
        histogram is normalised to fractional weights (sum = 1) so different
        snapshots or weight arrays are directly comparable.

        Parameters
        ----------
        weights : (N³,) ndarray, optional
            Per-cell weight array.  Defaults to uniform cell volume (vol_cell),
            giving a volume-weighted PDF.  Pass e.g. a line-luminosity array to
            get an emission-weighted phase diagram.
        bins : int or (int, int)
            Number of bins along the ne and Te axes.  Default 128.
        ne_range : (float, float), optional
            Linear limits for the ne axis [cm⁻³], e.g. (0.1, 1e4).  Defaults to
            data range.
        Te_range : (float, float), optional
            Linear limits for the Te axis [K], e.g. (5e3, 1e5).  Defaults to
            data range.
        mask_ionized : bool
            Restrict to ionized cells (x_HII > 0.2).  Default True.
        cmap : str, optional
            Colormap.  Defaults to 'viridis'.
        norm_log : bool
            Use logarithmic color normalization.  Default True.
        cbar_label : str, optional
            Colorbar label.  Auto-set based on weight type if omitted.
        fig, ax : optional
            Existing figure/axes to draw into.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ax  : matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.gridspec import GridSpec

        mask = self.ionized_mask() if mask_ionized else np.ones(self.N ** 3, dtype=bool)

        ne_all = self.ne
        Te_all = self.Te

        if weights is None:
            w_all = np.full(ne_all.shape, self.vol_cell)
            _cbar_label = cbar_label or 'Volume fraction'
        else:
            w_all = np.asarray(weights, dtype=float)
            _cbar_label = cbar_label or 'Weighted fraction'

        ne_sel = ne_all[mask]
        Te_sel = Te_all[mask]
        w_sel  = w_all[mask]

        valid  = (ne_sel > 0) & (Te_sel > 0) & np.isfinite(w_sel)
        ne_v   = ne_sel[valid]
        Te_v   = Te_sel[valid]
        w      = w_sel[valid]

        ne_lim = ne_range if ne_range is not None else (float(ne_v.min()), float(ne_v.max()))
        Te_lim = Te_range if Te_range is not None else (float(Te_v.min()), float(Te_v.max()))

        nb     = (bins, bins) if isinstance(bins, int) else tuple(bins)
        xe     = np.geomspace(ne_lim[0], ne_lim[1], nb[0] + 1)
        ye     = np.geomspace(Te_lim[0], Te_lim[1], nb[1] + 1)
        H, _, _ = np.histogram2d(ne_v, Te_v, bins=[xe, ye], weights=w)
        total  = H.sum()
        if total > 0:
            H /= total

        # ── Layout ───────────────────────────────────────────────────────
        fresh = (fig is None or ax is None)
        if fresh:
            fig = plt.figure(figsize=(7, 7))
            gs  = GridSpec(2, 2, figure=fig,
                           width_ratios=[1, 0.3], height_ratios=[0.3, 1],
                           hspace=0.0, wspace=0.0)
            ax      = fig.add_subplot(gs[1, 0])
            ax_top  = fig.add_subplot(gs[0, 0], sharex=ax)
            ax_side = fig.add_subplot(gs[1, 1], sharey=ax)
            fig.add_subplot(gs[0, 1]).set_visible(False)
        else:
            ax_top = ax_side = None

        # ── 2D histogram ─────────────────────────────────────────────────
        _cmap    = cmap or 'viridis'
        H_masked = np.ma.masked_where(H == 0, H)
        pos      = H_masked.compressed()
        if norm_log and len(pos):
            vmax = float(pos.max())
            vmin = vmax * 1e-3
            norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = None

        im = ax.pcolormesh(xe, ye, H_masked.T, cmap=_cmap, norm=norm)

        if cbar:
            cax = ax.inset_axes([0.57, 0.90, 0.40, 0.04])
            cb  = fig.colorbar(im, cax=cax, orientation='horizontal')
            cb.set_label(_cbar_label, fontsize=14, labelpad=2)
            cb.ax.tick_params(labelsize=12, top=True, bottom=False,
                              labeltop=True, labelbottom=False)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$n_e\;[\mathrm{cm}^{-3}]$', fontsize=18)
        ax.set_ylabel(r'$T_e\;[\mathrm{K}]$', fontsize=18)

        # ── Marginal 1D PDFs ──────────────────────────────────────────────
        if fresh:
            ne_centers = np.sqrt(xe[:-1] * xe[1:])   # geometric bin centres
            Te_centers = np.sqrt(ye[:-1] * ye[1:])
            ne_pdf = H.sum(axis=1)   # marginal over Te
            Te_pdf = H.sum(axis=0)   # marginal over ne

            ax_top.plot(ne_centers, ne_pdf, color='k', lw=1)
            ax_top.set_xscale('log')
            ax_top.set_yscale('log')
            ax_top.set_ylim(bottom=ne_pdf.max() * 1e-5)
            ax_top.tick_params(labelbottom=False, labelsize=14)
            ax_top.set_ylabel('PDF', fontsize=14)
            ax_top.set_xlim(ne_lim)

            ax_side.plot(Te_pdf, Te_centers, color='k', lw=1)
            ax_side.set_xscale('log')
            ax_side.set_yscale('log')
            ax_side.set_xlim(left=Te_pdf.max() * 1e-3)
            ax_side.tick_params(labelleft=False, labelsize=14)
            ax_side.set_xlabel('PDF', fontsize=14)
            ax_side.set_ylim(Te_lim)

        return fig, ax

    # ── Resource management ───────────────────────────────────────────────────

    def close(self):
        """Close the underlying HDF5 file handles."""
        self._states.close()
        self._colt.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self) -> str:
        return (f"COLTSnapshot(snapnum={self.snapnum}, "
                f"Zgas_Zsun={self.Zgas_Zsun}, N={self.N}, "
                f"thermal_eq={self.thermal_eq})")

    # ── Internal ──────────────────────────────────────────────────────────────

    def _load(self, key: str, loader) -> object:
        if key not in self._cache:
            self._cache[key] = loader()
        return self._cache[key]
