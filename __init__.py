"""
nebular_em — pyneb-based nebular emission line luminosity and diagnostics.

Compute per-cell or integrated line luminosities from simulation grid data
(temperature, electron density, ion number densities), then derive physical
conditions (Te, ne) from line ratios using pyneb.  Luminosities from external
codes (e.g. COLT) can be fed directly into the diagnostics layer.

Quick start
-----------
>>> from nebular_em import GridEmission, Diagnostics, LINE_LIBRARY
>>> grid = GridEmission(Te, ne,
...                     ion_densities={'OIII': n_OIII, 'OII': n_OII,
...                                    'HII': n_HII,  'SII': n_SII},
...                     vol_cell=vol_cell)
>>> L = grid.integrated(['OIII_5007', 'OIII_4363', 'OII_3726', 'OII_3729', 'Ha', 'Hb'])
>>> diag = Diagnostics(L)
>>> Te, ne = diag.iterate_Te_ne(('OIII_4363', 'OIII_5007'), ('OII_3726', 'OII_3729'))
"""

from .emission    import GridEmission, LineSpec
from .diagnostics import Diagnostics
from .line_library import LINE_LIBRARY, LINE_LABELS, LINE_ALIASES
from .colt        import COLTSnapshot, ION_MAP, SOLAR_MASS_FRAC, solar_number_abundance

__all__    = ['GridEmission', 'LineSpec', 'Diagnostics', 'LINE_LIBRARY', 'LINE_LABELS', 'LINE_ALIASES',
              'COLTSnapshot', 'ION_MAP', 'SOLAR_MASS_FRAC', 'solar_number_abundance']
__version__ = '0.1.0'
