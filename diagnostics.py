"""
diagnostics.py — derive electron temperature and density from line ratios.

The Diagnostics class wraps pyneb's getTemDen() for named line pairs from
LINE_LIBRARY.  Luminosities may come from GridEmission.integrated() or from
an external code such as COLT — the interface is the same.
"""

from __future__ import annotations
import numpy as np

# OII V1 recombination multiplet keys used for O²⁺/H via recombination lines.
# These 6 lines match the set validated in the notebook; 4674 and 4696 omitted
# as they have weaker / less reliable atomic data.
_OII_REC_KEYS: list[str] = [
    'OII_rec_4639', 'OII_rec_4642', 'OII_rec_4649',
    'OII_rec_4651', 'OII_rec_4662', 'OII_rec_4676',
]


def _safe(v: object) -> float:
    """Return float(v) if finite and positive, else nan.  Guards all pyneb getTemDen calls."""
    try:
        v = float(v)
        return v if (np.isfinite(v) and v > 0) else float('nan')
    except Exception:
        return float('nan')


class Diagnostics:
    """
    Derive Te and ne from integrated emission line ratios using pyneb.

    Parameters
    ----------
    luminosities : dict[str, float]
        Total line luminosities [erg/s], keyed by LINE_LIBRARY names (or any
        consistent naming scheme matching what you pass to the diagnostic methods).
        Extra keys are silently ignored.

    Examples
    --------
    >>> diag = Diagnostics({'OII_3726': 1e38, 'OII_3729': 2e38,
    ...                     'OIII_4363': 3e36, 'OIII_5007': 5e39})
    >>> ne = diag.electron_density('OII_3726', 'OII_3729', Te=1e4)
    >>> Te = diag.electron_temperature('OIII_4363', 'OIII_5007', ne=ne)
    >>> Te, ne = diag.iterate_Te_ne(('OIII_4363', 'OIII_5007'),
    ...                              ('OII_3726',  'OII_3729'))
    """

    def __init__(self, luminosities: dict[str, float]):
        self.L = dict(luminosities)

    # ------------------------------------------------------------------
    # Public diagnostics
    # ------------------------------------------------------------------

    def electron_density(self, line1: str, line2: str, Te: float = 1e4) -> float:
        """
        Derive ne from L[line1] / L[line2] at an assumed electron temperature.

        Both lines must belong to the same pyneb atom (same ion species).

        Parameters
        ----------
        line1, line2 : LINE_LIBRARY keys for the numerator and denominator lines.
        Te           : assumed electron temperature [K].

        Returns
        -------
        float — electron density [cm⁻³]
        """
        s1, s2 = self._spec_pair(line1, line2)
        return _safe(s1.atom.getTemDen(
            self.L[line1] / self.L[line2],
            tem=Te,
            to_eval=self._to_eval(s1, s2),
        ))

    def electron_temperature(self, line1: str, line2: str, ne: float = 100.) -> float:
        """
        Derive Te from L[line1] / L[line2] at an assumed electron density.

        Both lines must belong to the same pyneb atom.

        Parameters
        ----------
        line1, line2 : LINE_LIBRARY keys for the numerator and denominator lines.
        ne           : assumed electron density [cm⁻³].

        Returns
        -------
        float — electron temperature [K], or nan if pyneb cannot converge.
        """
        s1, s2 = self._spec_pair(line1, line2)
        return _safe(s1.atom.getTemDen(
            self.L[line1] / self.L[line2],
            den=ne,
            to_eval=self._to_eval(s1, s2),
        ))

    def iterate_Te_ne(
        self,
        Te_lines: tuple[str, str],
        ne_lines: tuple[str, str],
        Te0: float = 1e4,
        ne0: float = 100.,
        n_iter: int = 5,
    ) -> tuple[float, float]:
        """
        Self-consistently solve for Te and ne by iterating two diagnostics.

        Alternates between the temperature diagnostic (holding ne fixed) and the
        density diagnostic (holding Te fixed) until convergence.  3-5 iterations
        are typically sufficient.

        Parameters
        ----------
        Te_lines : (line1, line2) for the temperature ratio, e.g. ('OIII_4363', 'OIII_5007')
        ne_lines : (line1, line2) for the density  ratio, e.g. ('OII_3726',  'OII_3729')
        Te0, ne0 : initial guesses [K, cm⁻³]
        n_iter   : number of alternating iterations

        Returns
        -------
        (Te, ne) : converged electron temperature [K] and density [cm⁻³]
        """
        Te, ne = Te0, ne0
        for _ in range(n_iter):
            Te = self.electron_temperature(*Te_lines, ne=ne)
            if not np.isfinite(Te):
                return float('nan'), float('nan')
            ne = self.electron_density(*ne_lines, Te=Te)
            if not np.isfinite(ne):
                return float('nan'), float('nan')
        return Te, ne

    # ------------------------------------------------------------------
    # Per-cell diagnostics (require a GridEmission object)
    # ------------------------------------------------------------------

    @staticmethod
    def temperature_variance(
        grid: 'GridEmission',
        weight_line: str,
        mask=None,
    ) -> tuple[float, float]:
        """
        Emission-measure–weighted mean temperature T0 and Peimbert's t² for a
        line-emitting region.

            w_i = n_{e,i} × n_{ion,i}          (emission-measure weight)
            T0  = Σ_i w_i T_i / Σ_i w_i
            t²  = Σ_i w_i (T_i − T0)² / (T0² Σ_i w_i)

        The ion species is taken from the weight_line's ion_key in LINE_LIBRARY
        (e.g. 'OIII_4363' → n_OIII, 'NII_5755' → n_NII).  Cell volume cancels
        for a uniform grid and is not included in the weight.

        Parameters
        ----------
        grid        : GridEmission
            Pre-built GridEmission (e.g. from COLTSnapshot.to_grid_emission).
        weight_line : str
            LINE_LIBRARY key identifying the ion zone, e.g. 'OIII_4363'.
        mask        : array-like of bool, optional
            Cell selection mask (e.g. snap.ionized_mask()).

        Returns
        -------
        T0 : float — emission-measure–weighted mean electron temperature [K]
        t2 : float — Peimbert temperature variance (dimensionless)
        """
        from .line_library import LINE_LIBRARY
        ion_key = LINE_LIBRARY[weight_line].ion_key   # e.g. 'OIII' or 'NII'

        Te    = grid.Te
        ne    = grid.ne
        n_ion = grid.ions[ion_key]
        w     = ne * n_ion   # emission-measure weight; uniform V cancels

        if mask is not None:
            mask = np.asarray(mask)
            idx  = np.where(mask)[0] if mask.dtype == bool else mask.astype(int)
            Te    = Te[idx]
            w     = w[idx]

        pos = w > 0
        if pos.sum() < 2:
            return float('nan'), float('nan')

        T0 = float(np.average(Te[pos], weights=w[pos]))
        t2 = float(np.average((Te[pos] - T0) ** 2, weights=w[pos]) / T0 ** 2)
        return T0, t2

    def get_Oabundance(
        self,
        Te: float,
        ne: float,
        type: str = 'O',
        method: str = 'cel',
    ) -> float:
        """
        Derive oxygen ionic or total abundance relative to hydrogen.

        Uses the direct method (Te, ne known) with the two-ionization-zone
        approximation O/H ≈ O⁺/H + O²⁺/H.

        Parameters
        ----------
        Te : float
            Electron temperature [K], e.g. from iterate_Te_ne.
        ne : float
            Electron density [cm⁻³].
        type : {'O', 'OII', 'OIII'}
            'OII'  — O⁺/H  (always CEL: [O II] 3726+3729 / Hβ)
            'OIII' — O²⁺/H (CEL: [O III] 4959+5007 / Hβ;
                             REC: OII V1 multiplet ~4639–4676 Å / Hβ)
            'O'    — total (O⁺ + O²⁺) / H
        method : {'cel', 'rec'}
            'cel' — collisionally excited lines.
            'rec' — O²⁺/H from OII V1 recombination multiplet;
                    O⁺/H still uses CEL (no suitable recombination lines).

        Returns
        -------
        float — linear abundance n(O)/n(H), or nan if required luminosities
                are missing or non-positive.  Convert to log-scale with
                12 + np.log10(result).
        """
        from .line_library import LINE_LIBRARY

        if type not in ('O', 'OII', 'OIII'):
            raise ValueError(f"type must be 'O', 'OII', or 'OIII'; got '{type!r}'")
        if method not in ('cel', 'rec'):
            raise ValueError(f"method must be 'cel' or 'rec'; got '{method!r}'")

        _o3  = LINE_LIBRARY['OIII_5007'].atom
        _o2  = LINE_LIBRARY['OII_3726'].atom
        _o2r = LINE_LIBRARY['OII_rec_4639'].atom

        L  = self.L
        Hb = L.get('Hb', 0.0)
        if not (np.isfinite(Hb) and Hb > 0):
            return float('nan')

        def _OII_cel() -> float:
            L3726 = L.get('OII_3726', 0.0)
            L3729 = L.get('OII_3729', 0.0)
            if not (L3726 + L3729 > 0):
                return float('nan')
            return _safe(_o2.getIonAbundance(
                (L3726 + L3729) / Hb * 100,
                Te, ne, to_eval='L(3726)+L(3729)',
            ))

        def _OIII_cel() -> float:
            L4959 = L.get('OIII_4959', 0.0)
            L5007 = L.get('OIII_5007', 0.0)
            if not (L4959 + L5007 > 0):
                return float('nan')
            return _safe(_o3.getIonAbundance(
                (L4959 + L5007) / Hb * 100,
                Te, ne, to_eval='L(4959)+L(5007)',
            ))

        def _OIII_rec() -> float:
            avail = [k for k in _OII_REC_KEYS if L.get(k, 0.0) > 0]
            if not avail:
                return float('nan')
            lum_rec = sum(L[k] for k in avail)
            to_eval = ' + '.join(f'L({LINE_LIBRARY[k].wav})' for k in avail)
            return _safe(_o2r.getIonAbundance(
                lum_rec / Hb * 100,
                Te, ne, to_eval=to_eval,
            ))

        if type == 'OII':
            return _OII_cel()
        if type == 'OIII':
            return _OIII_cel() if method == 'cel' else _OIII_rec()

        # type == 'O': total
        oii  = _OII_cel()
        oiii = _OIII_cel() if method == 'cel' else _OIII_rec()
        if not (np.isfinite(oii) and np.isfinite(oiii)):
            return float('nan')
        return oii + oiii

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _spec_pair(self, line1: str, line2: str):
        """Look up two LineSpecs and verify they share the same pyneb atom."""
        from .line_library import LINE_LIBRARY
        for k in (line1, line2):
            if k not in LINE_LIBRARY:
                raise KeyError(
                    f"'{k}' not in LINE_LIBRARY.  For custom lines use pyneb directly."
                )
        s1, s2 = LINE_LIBRARY[line1], LINE_LIBRARY[line2]
        if s1.atom is not s2.atom:
            raise ValueError(
                f"'{line1}' and '{line2}' belong to different pyneb atoms "
                f"and cannot be combined in a single getTemDen call."
            )
        return s1, s2

    @staticmethod
    def _to_eval(s1, s2) -> str:
        """Build the pyneb to_eval string for a ratio of two transitions."""
        return f"I({s1.lev_i}, {s1.lev_j}) / I({s2.lev_i}, {s2.lev_j})"
