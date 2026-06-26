# nebular_em

A Python library for computing nebular emission line luminosities from simulation grid data and deriving physical conditions (T_e, n_e) from line ratios using [PyNeb](https://www.iaa.csic.es/~wetal/pyneb/).

**Why this library?**
- Compute luminosities for lines your radiative transfer code does not produce (e.g. UV OIII], CIII], NIV]).
- Cross-check code-output luminosities with an independent pyneb-based calculation.
- Perform observer-style diagnostics (T_e, n_e) from any set of integrated line luminosities — whether computed here or taken directly from your RT code.

---

## Dependencies

- Python ≥ 3.10
- [pyneb](https://pypi.org/project/PyNeb/) ≥ 1.1
- numpy, scipy, h5py

---

## Installation

```bash
# from the repo root
pip install -e nebular_em/
# or simply add the repo root to PYTHONPATH:
export PYTHONPATH=/path/to/nebular-emission:$PYTHONPATH
```

---

## Quick start

### With COLT/FLASH snapshots

`COLTSnapshot` handles all the file I/O and ion-density conversion, then hands off to the core classes:

```python
from nebular_em import COLTSnapshot, Diagnostics

snap_dir = '/path/to/colt_dir'   # directory with states-teq_NNNN.hdf5 + colt_NNNN.hdf5

with COLTSnapshot(snap_dir, snapnum=150, Zgas_Zsun=0.1) as snap:
    # Build a GridEmission from the snapshot
    grid = snap.to_grid_emission(['OIII', 'OII', 'HII', 'SII'])

    # Compute integrated luminosities over the ionized region
    lines = ['OIII_5007', 'OIII_4363', 'OII_3726', 'OII_3729', 'Ha', 'Hb']
    L = grid.integrated(lines, mask=snap.ionized_mask())

    # Observer-style diagnostics
    Te, ne = Diagnostics(L).iterate_Te_ne(
        Te_lines=('OIII_4363', 'OIII_5007'),
        ne_lines=('OII_3726',  'OII_3729'),
    )
    print(f"T_e = {Te:.0f} K,  n_e = {ne:.1f} cm⁻³")
```

### With arbitrary simulation data

```python
from nebular_em import GridEmission, Diagnostics, LINE_LIBRARY

# Te, n_e, n_OIII, ... are (N,) arrays; vol_cell is a float or (N,) array
grid = GridEmission(
    Te=Te,
    ne=n_e,
    ion_densities={'OIII': n_OIII, 'OII': n_OII, 'SII': n_SII, 'HII': n_HII},
    vol_cell=vol_cell,
)

lines    = ['OIII_5007', 'OIII_4363', 'OII_3726', 'OII_3729', 'SII_6716', 'SII_6731', 'Ha', 'Hb']
cell_L   = grid.cell_luminosities(lines, mask=ionized)   # dict of (N,) arrays [erg/s]
total_L  = grid.integrated(lines, mask=ionized)          # dict of floats [erg/s]

diag = Diagnostics(total_L)    # also accepts COLT-output luminosity dicts
Te, ne = diag.iterate_Te_ne(('OIII_4363', 'OIII_5007'), ('OII_3726', 'OII_3729'))
```

---

## API reference

### `COLTSnapshot`

```python
COLTSnapshot(snap_dir, snapnum, Zgas_Zsun,
             N=512, X_H=0.76, thermal_eq=True)
```

Reads `states-teq_NNNN.hdf5` and `colt_NNNN.hdf5` from `snap_dir`.  All arrays are loaded lazily and cached on first access.  Use as a context manager (`with COLTSnapshot(...) as snap`) to ensure HDF5 file handles are closed.

| Parameter | Description |
|---|---|
| `snap_dir` | Directory containing the snapshot HDF5 files |
| `snapnum` | Snapshot number (integer), used to build filenames |
| `Zgas_Zsun` | Uniform gas metallicity in solar units — must match the COLT run |
| `N` | Grid cells per side (default 512) |
| `X_H` | Hydrogen mass fraction (default 0.76) |
| `thermal_eq` | If `True` (default) load `states-teq_NNNN.hdf5`; if `False` load `states_NNNN.hdf5` |

#### Properties

| Property | Description |
|---|---|
| `Te` | Equilibrium electron temperature [K], shape (N³,) |
| `nH` | Hydrogen nucleon density [cm⁻³] = nHI + nHII + 2 nH₂ |
| `ne` | Electron density [cm⁻³] from post-COLT ionization state |
| `x_HII` | HII ionization fraction = n_HII / n_H_nuclei |
| `vol_cell` | Cell volume [cm³] (uniform cubic grid) |
| `bbox` | Bounding box [cm], shape (2, 3) |

#### Methods

**`ion_densities(species) → dict[str, ndarray]`**

Compute number densities [cm⁻³] for a list of ion species using solar abundances (Grevesse et al. 2010) scaled by `Zgas_Zsun`.  Results are cached.  See `ION_MAP` for all supported keys.

```python
ions = snap.ion_densities(['OIII', 'OII', 'HII', 'SII', 'CIII', 'NIV', 'NeIII'])
```

**`ionized_mask(threshold=0.2) → ndarray`**

Boolean mask of cells with `x_HII > threshold`.

**`to_grid_emission(species, **kwargs) → GridEmission`**

Convenience wrapper: calls `ion_densities(species)` and constructs a `GridEmission`.  `**kwargs` are forwarded to `GridEmission` (e.g. `Te_range`, `ne_range`).

---

### `GridEmission`

```python
GridEmission(Te, ne, ion_densities, vol_cell,
             Te_range=(3e3, 1e5), ne_range=(1e-2, 1e7),
             n_Te=80, n_ne=80)
```

| Parameter | Type | Description |
|---|---|---|
| `Te` | `(N,) array` | Electron temperature [K] |
| `ne` | `(N,) array` | Electron density [cm⁻³] |
| `ion_densities` | `dict[str, (N,) array]` | Ion number densities [cm⁻³], e.g. `{'OIII': n_OIII}` |
| `vol_cell` | `float` or `(N,) array` | Cell volume [cm³] |
| `Te_range` | `(float, float)` | Temperature bounds of the interpolation table [K] |
| `ne_range` | `(float, float)` | Density bounds of the interpolation table [cm⁻³] |
| `n_Te`, `n_ne` | `int` | Grid resolution (log-spaced) along each axis |

Emissivity tables are built lazily on first request per line and cached for the lifetime of the object.

#### Methods

**`cell_luminosities(lines, mask=None) → dict[str, ndarray]`**

Per-cell luminosity [erg/s].  `lines` is a list of `LINE_LIBRARY` keys or a custom `{name: LineSpec}` dict.  Cells outside `mask` are zero.

**`integrated(lines, mask=None) → dict[str, float]`**

Sum of `cell_luminosities` over all (masked) cells.

---

### `Diagnostics`

```python
Diagnostics(luminosities: dict[str, float])
```

Accepts any `{line_name: L_total [erg/s]}` dict — from `GridEmission.integrated()` or directly from an external code.

#### Methods

**`electron_density(line1, line2, Te=1e4) → float`**

Derive n_e [cm⁻³] from L[line1]/L[line2] at the assumed T_e.
Both lines must be from the same ion (same pyneb atom).

**`electron_temperature(line1, line2, ne=100.) → float`**

Derive T_e [K] from L[line1]/L[line2] at the assumed n_e.

**`iterate_Te_ne(Te_lines, ne_lines, Te0=1e4, ne0=100., n_iter=5) → (float, float)`**

Iterates the T_e and n_e diagnostics to self-consistency.  Returns `(Te, ne)`.

---

### `LineSpec`

```python
LineSpec(atom, lev_i, lev_j, ion_key, wav, rec=False)
```

Frozen dataclass describing one emission line.  Used to extend `LINE_LIBRARY` with custom lines (see [Adding a custom line](#adding-a-custom-line)).

---

## Built-in lines (`LINE_LIBRARY`)

| Key | Ion | λ (Å) | Levels | Diagnostic use |
|---|---|---|---|---|
| `OIII_5007` | O²⁺ | 5006.8 | 4→3 | BPT, strong-line O/H |
| `OIII_4363` | O²⁺ | 4363.2 | 5→4 | T_e diagnostic |
| `OIII_88` | O²⁺ | 883 564 (88.3 μm) | 2→1 | n_e diagnostic (FIR) |
| `OIII_52` | O²⁺ | 518 145 (51.8 μm) | 3→2 | n_e diagnostic (FIR) |
| `OIII_1661` | O²⁺ | 1660.8 | 6→2 | n_e diagnostic (UV, AK99) |
| `OIII_1666` | O²⁺ | 1666.1 | 6→3 | n_e diagnostic (UV, AK99) |
| `OII_3726` | O⁺ | 3726.0 | 3→1 | n_e diagnostic |
| `OII_3729` | O⁺ | 3728.8 | 2→1 | n_e diagnostic |
| `SII_6716` | S⁺ | 6716.4 | 3→1 | n_e diagnostic |
| `SII_6731` | S⁺ | 6730.8 | 2→1 | n_e diagnostic |
| `CIII_1907` | C²⁺ | 1906.7 | 4→1 | n_e diagnostic (UV) |
| `CIII_1909` | C²⁺ | 1908.7 | 3→1 | n_e diagnostic (UV) |
| `NIV_1483` | N³⁺ | 1483.3 | 4→1 | n_e diagnostic (UV) |
| `NIV_1486` | N³⁺ | 1486.5 | 3→1 | n_e diagnostic (UV) |
| `NII_6583` | N⁺ | 6583.4 | 4→3 | BPT |
| `NII_5755` | N⁺ | 5754.6 | 5→4 | T_e diagnostic |
| `Ha` | H⁺ | 6562.8 | 3→2 | SFR tracer, BPT denominator |
| `Hb` | H⁺ | 4861.3 | 4→2 | SFR tracer, BPT denominator |
| `Hg` | H⁺ | 4340.5 | 5→2 | Balmer decrement |

The UV OIII lines (1661, 1666 Å) use Aggarwal & Keenan (1999) collision strengths via a separate pyneb atom instance.

---

## Ion species supported by `COLTSnapshot` (`ION_MAP`)

All ionization stages tracked in the COLT states file are available:

| Group | Keys |
|---|---|
| Hydrogen | `HI`, `HII` |
| Helium | `HeI`, `HeII` |
| Carbon | `CI`, `CII`, `CIII`, `CIV` |
| Nitrogen | `NI`, `NII`, `NIII`, `NIV`, `NV` |
| Oxygen | `OI`, `OII`, `OIII`, `OIV` |
| Neon | `NeI`, `NeII`, `NeIII`, `NeIV` |
| Magnesium | `MgI`, `MgII`, `MgIII` |
| Silicon | `SiI`, `SiII`, `SiIII`, `SiIV` |
| Sulfur | `SI`, `SII`, `SIII`, `SIV`, `SV` |
| Iron | `FeI`, `FeII`, `FeIII`, `FeIV`, `FeV`, `FeVI` |

Number densities are computed as:
- **H**: `nH × x`
- **He**: `nH × x × X_He / (4 X_H)` where X_He = 1 − X_H
- **Metals**: `nH × x × (n_X/n_H)_solar × Zgas_Zsun`

with solar abundances from Grevesse et al. (2010).

---

## Adding a custom line

Use `LineSpec` directly and pass a `{name: LineSpec}` dict to `cell_luminosities` / `integrated`:

```python
import pyneb as pn
from nebular_em import GridEmission, LineSpec

# Example: [Ne III] 3869 Å (lev_i=3, lev_j=1)
ne3 = pn.Atom('Ne', 3)
custom_lines = {
    'NeIII_3869': LineSpec(atom=ne3, lev_i=3, lev_j=1, ion_key='NeIII', wav=3869.),
}
# ion_densities must include 'NeIII'
grid = GridEmission(Te, ne, ion_densities={'NeIII': n_NeIII}, vol_cell=vol_cell)
cell_L = grid.cell_luminosities(custom_lines)
```

With COLT data, `ion_densities` is already provided by `COLTSnapshot`:
```python
with COLTSnapshot(snap_dir, 150, Zgas_Zsun=0.1) as snap:
    grid = snap.to_grid_emission(['NeIII'])
    cell_L = grid.cell_luminosities(custom_lines, mask=snap.ionized_mask())
```

To use a custom line in `Diagnostics`, call pyneb directly:

```python
ratio = total_L['NeIII_3869'] / total_L['NeIII_3967']
ne3.getTemDen(ratio, tem=1e4, to_eval="I(3, 1) / I(2, 1)")
```

---

## Notes

### Emissivity interpolation

- Tables are log-spaced in both T_e and n_e and interpolated bilinearly in log–log space.
- `fill_value=None` means the interpolator extrapolates outside the grid bounds; set `Te_range`/`ne_range` to cover your data.
- Some pyneb recombination tables (e.g. Storey & Hummer 1995 for H I) only cover T_e ≲ 30 000 K.  Cells outside this range are filled by nearest-neighbour extrapolation from the table boundary; errors are small for lines whose emissivity is a smooth power law in T_e.
- Processing is chunked at 5 × 10⁶ cells to limit peak memory.
- Emissivity tables are cached per `LineSpec` instance — instantiate `GridEmission` once and call `cell_luminosities` / `integrated` multiple times at no extra cost.

### Computing ion densities without `COLTSnapshot`

If you are not using COLT, convert simulation ion fractions manually:

```python
PROTONMASS = 1.6726e-24   # g
X_H = 0.76

nH     = rho * X_H / PROTONMASS             # hydrogen nucleon density [cm⁻³]
O_solar = 0.00573805 * 1.008 / (X_H * 15.999)   # n_O/n_H at solar Z (Grevesse+2010)
n_OIII  = x_OIII * nH * O_solar * (Z_gas / Z_sun)
```

The `solar_number_abundance(element, X_H)` helper from `nebular_em` performs this calculation for any element in `SOLAR_MASS_FRAC`.
