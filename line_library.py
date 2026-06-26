"""
LINE_LIBRARY — pre-defined emission line specifications using shared pyneb atom instances.

Atom instances are created once at module import.  The UV OIII lines (1661, 1666 Å)
use the Aggarwal & Keenan (1999) collision strengths via a separate atom instance;
pyneb's global data file pointer is restored to defaults afterward.
"""

import pyneb as pn
from .emission import LineSpec

# ---------------------------------------------------------------------------
# Atom instances (shared across all LineSpecs referencing the same ion)
# ---------------------------------------------------------------------------
_o3 = pn.Atom('O', 3)
_o2 = pn.Atom('O', 2)
_s2 = pn.Atom('S', 2)
_n2 = pn.Atom('N', 2)
_c3 = pn.Atom('C', 3)
_c4 = pn.Atom('C', 4)
_n4 = pn.Atom('N', 4)
_H1    = pn.RecAtom('H', 1, case='B')
_o2_rec = pn.RecAtom('O', 2)   # OII recombination: O²⁺ + e⁻ → O⁺ + hν; emissivity ∝ n_OIII × n_e

# UV OIII: Aggarwal & Keenan (1999) collision strengths needed for 1661/1666 Å
pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
_o3_uv = pn.Atom('O', 3)
try:
    pn.atomicData.resetDataFileDict()   # restore defaults so we don't pollute user's pyneb state
except AttributeError:
    pass

# ---------------------------------------------------------------------------
# LINE_LIBRARY
# Keys follow the pattern  ION_WAVELENGTH  where wavelength is in Å
# (FIR lines quoted in μm for clarity; stored wavelength is still in Å).
# ---------------------------------------------------------------------------
LINE_LIBRARY: dict[str, LineSpec] = {

    # ── [O III] ─────────────────────────────────────────────────────────────
    'OIII_5007': LineSpec(_o3,    lev_i=4, lev_j=3, ion_key='OIII', wav=5006.8),
    'OIII_4959': LineSpec(_o3,  lev_i=4, lev_j=2, ion_key='OIII', wav=4958.9),  # 4959
    'OIII_4363': LineSpec(_o3,    lev_i=5, lev_j=4, ion_key='OIII', wav=4363.2),
    'OIII_88'  : LineSpec(_o3,    lev_i=2, lev_j=1, ion_key='OIII', wav=883564.),   # 88.3 μm
    'OIII_52'  : LineSpec(_o3,    lev_i=3, lev_j=2, ion_key='OIII', wav=518145.),   # 51.8 μm
    'OIII_1661': LineSpec(_o3_uv, lev_i=6, lev_j=2, ion_key='OIII', wav=1660.8),   # UV, AK99
    'OIII_1666': LineSpec(_o3_uv, lev_i=6, lev_j=3, ion_key='OIII', wav=1666.1),   # UV, AK99

    # ── [O II] ──────────────────────────────────────────────────────────────
    'OII_3726'  : LineSpec(_o2,   lev_i=3, lev_j=1, ion_key='OII',  wav=3726.0),
    'OII_3729'  : LineSpec(_o2,   lev_i=2, lev_j=1, ion_key='OII',  wav=3728.8),
    'OII_7320'  : LineSpec(_o2, lev_i=5, lev_j=2, ion_key='OII',  wav=7318.9),
    'OII_7330'  : LineSpec(_o2, lev_i=5, lev_j=3, ion_key='OII',  wav=7329.7),
    # ── [S II] ──────────────────────────────────────────────────────────────
    'SII_6716'  : LineSpec(_s2,   lev_i=3, lev_j=1, ion_key='SII',  wav=6716.4),
    'SII_6731'  : LineSpec(_s2,   lev_i=2, lev_j=1, ion_key='SII',  wav=6730.8),

    # ── C III] / [C III] ────────────────────────────────────────────────────
    'CIII_1907' : LineSpec(_c3,   lev_i=4, lev_j=1, ion_key='CIII', wav=1906.7),
    'CIII_1909' : LineSpec(_c3,   lev_i=3, lev_j=1, ion_key='CIII', wav=1908.7),
    'CIV_1548'  : LineSpec(_c4, lev_i=3, lev_j=1, ion_key='CIV',  wav=1548.2),  # resonance doublet
    'CIV_1551'  : LineSpec(_c4, lev_i=2, lev_j=1, ion_key='CIV',  wav=1550.8),

    # ── N IV] / [N IV] ──────────────────────────────────────────────────────
    'NIV_1483'  : LineSpec(_n4,   lev_i=4, lev_j=1, ion_key='NIV',  wav=1483.3),
    'NIV_1486'  : LineSpec(_n4,   lev_i=3, lev_j=1, ion_key='NIV',  wav=1486.5),

    # ── [N II] ──────────────────────────────────────────────────────────────
    'NII_6583'  : LineSpec(_n2,   lev_i=4, lev_j=3, ion_key='NII',  wav=6583.4),
    'NII_6548'  : LineSpec(_n2, lev_i=4, lev_j=2, ion_key='NII',  wav=6548.1),
    'NII_5755'  : LineSpec(_n2,   lev_i=5, lev_j=4, ion_key='NII',  wav=5754.6),

    # ── O II recombination (V1 multiplet, ~4639–4696 Å) ─────────────────────
    # RecAtom('O', 2): O²⁺ + e⁻ → O⁺ + hν, emissivity ∝ n_OIII × n_e
    'OII_rec_4639': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4638.86, rec=True),
    'OII_rec_4642': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4641.81, rec=True),
    'OII_rec_4649': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4649.13, rec=True),
    'OII_rec_4651': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4650.84, rec=True),
    'OII_rec_4662': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4661.63, rec=True),
    'OII_rec_4674': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4673.73, rec=True),
    'OII_rec_4676': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4676.23, rec=True),
    'OII_rec_4696': LineSpec(_o2_rec, lev_i=None, lev_j=None, ion_key='OIII', wav=4696.35, rec=True),

    # ── H I recombination (case B) ───────────────────────────────────────────
    'Ha'        : LineSpec(_H1,   lev_i=3, lev_j=2, ion_key='HII',  wav=6562.8, rec=True),
    'Hb'        : LineSpec(_H1,   lev_i=4, lev_j=2, ion_key='HII',  wav=4861.3, rec=True),
    'Hg'        : LineSpec(_H1,   lev_i=5, lev_j=2, ion_key='HII',  wav=4340.5, rec=True),
}

# ---------------------------------------------------------------------------
# LINE_LABELS — LaTeX display labels for every LINE_LIBRARY key.
# Bracket convention: [X II] = forbidden, X II] = semi-forbidden, X II = allowed/recombination.
# Compatible with matplotlib mathtext (no usetex required).
# ---------------------------------------------------------------------------
LINE_LABELS: dict[str, str] = {
    # ── [O III] ──────────────────────────────────────────────────────────────
    'OIII_5007': r"$[\mathrm{O\,III}]\,\lambda5007$",
    'OIII_4959': r"$[\mathrm{O\,III}]\,\lambda4959$",
    'OIII_4363': r"$[\mathrm{O\,III}]\,\lambda4363$",
    'OIII_88'  : r"$[\mathrm{O\,III}]\,88\,\mu\mathrm{m}$",
    'OIII_52'  : r"$[\mathrm{O\,III}]\,52\,\mu\mathrm{m}$",
    'OIII_1661': r"$\mathrm{O\,III]}\,\lambda1661$",   # semi-forbidden
    'OIII_1666': r"$\mathrm{O\,III]}\,\lambda1666$",   # semi-forbidden
    # ── [O II] ───────────────────────────────────────────────────────────────
    'OII_3726'  : r"$[\mathrm{O\,II}]\,\lambda3726$",
    'OII_3729'  : r"$[\mathrm{O\,II}]\,\lambda3729$",
    'OII_7320'  : r"$[\mathrm{O\,II}]\,\lambda7319$",
    'OII_7330'  : r"$[\mathrm{O\,II}]\,\lambda7330$",
    # ── [S II] ───────────────────────────────────────────────────────────────
    'SII_6716'  : r"$[\mathrm{S\,II}]\,\lambda6716$",
    'SII_6731'  : r"$[\mathrm{S\,II}]\,\lambda6731$",
    # ── C III], C IV ─────────────────────────────────────────────────────────
    'CIII_1907' : r"$\mathrm{C\,III]}\,\lambda1907$",  # semi-forbidden
    'CIII_1909' : r"$\mathrm{C\,III]}\,\lambda1909$",  # semi-forbidden
    'CIV_1548'  : r"$\mathrm{C\,IV}\,\lambda1548$",
    'CIV_1551'  : r"$\mathrm{C\,IV}\,\lambda1551$",
    # ── N IV], [N II] ────────────────────────────────────────────────────────
    'NIV_1483'  : r"$\mathrm{N\,IV]}\,\lambda1483$",   # semi-forbidden
    'NIV_1486'  : r"$\mathrm{N\,IV]}\,\lambda1486$",   # semi-forbidden
    'NII_6583'  : r"$[\mathrm{N\,II}]\,\lambda6583$",
    'NII_6548'  : r"$[\mathrm{N\,II}]\,\lambda6548$",
    'NII_5755'  : r"$[\mathrm{N\,II}]\,\lambda5755$",
    # ── O II recombination (V1 multiplet) ────────────────────────────────────
    'OII_rec_4639': r"$\mathrm{O\,II}\,\lambda4639\,(\mathrm{rec})$",
    'OII_rec_4642': r"$\mathrm{O\,II}\,\lambda4642\,(\mathrm{rec})$",
    'OII_rec_4649': r"$\mathrm{O\,II}\,\lambda4649\,(\mathrm{rec})$",
    'OII_rec_4651': r"$\mathrm{O\,II}\,\lambda4651\,(\mathrm{rec})$",
    'OII_rec_4662': r"$\mathrm{O\,II}\,\lambda4662\,(\mathrm{rec})$",
    'OII_rec_4674': r"$\mathrm{O\,II}\,\lambda4674\,(\mathrm{rec})$",
    'OII_rec_4676': r"$\mathrm{O\,II}\,\lambda4676\,(\mathrm{rec})$",
    'OII_rec_4696': r"$\mathrm{O\,II}\,\lambda4696\,(\mathrm{rec})$",
    # ── H I recombination ────────────────────────────────────────────────────
    'Ha': r"H$\alpha$",
    'Hb': r"H$\beta$",
    'Hg': r"H$\gamma$",
    # ── Aliases (doublets / multiplet sums) ──────────────────────────────────
    'OII_3726_3729'  : r"$[\mathrm{O\,II}]\,\lambda\lambda3726,3729$",
    'OIII_4959_5007' : r"$[\mathrm{O\,III}]\,\lambda\lambda4959,5007$",
    'SII_6716_6731'  : r"$[\mathrm{S\,II}]\,\lambda\lambda6716,6731$",
    'NII_6548_6583'  : r"$[\mathrm{N\,II}]\,\lambda\lambda6548,6583$",
    'OII_7320_7330'  : r"$[\mathrm{O\,II}]\,\lambda\lambda7319,7330$",
    'CIII_1907_1909' : r"$\mathrm{C\,III]}\,\lambda\lambda1907,1909$",
    'CIV_1548_1551'  : r"$\mathrm{C\,IV}\,\lambda\lambda1548,1551$",
    'NIV_1483_1486'  : r"$\mathrm{N\,IV]}\,\lambda\lambda1483,1486$",
}

# ---------------------------------------------------------------------------
# LINE_ALIASES — map a summed-line key to its constituent LINE_LIBRARY keys.
# Add new entries here to extend the available sums in plot_emission_map.
# ---------------------------------------------------------------------------
LINE_ALIASES: dict[str, list[str]] = {
    'OII_3726_3729'  : ['OII_3726',  'OII_3729' ],
    'OIII_4959_5007' : ['OIII_4959', 'OIII_5007'],
    'SII_6716_6731'  : ['SII_6716',  'SII_6731' ],
    'NII_6548_6583'  : ['NII_6548',  'NII_6583' ],
    'OII_7320_7330'  : ['OII_7320',  'OII_7330' ],
    'CIII_1907_1909' : ['CIII_1907', 'CIII_1909'],
    'CIV_1548_1551'  : ['CIV_1548',  'CIV_1551' ],
    'NIV_1483_1486'  : ['NIV_1483',  'NIV_1486' ],
}
