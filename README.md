# Wind farm connectivity in the German EEZ – OceanTracker illustration

An example that uses OceanTracker with the SCHISM GermanBight hindcast
(`/hpcfreenas/hindcast/Europe/GermanBight`, January 2020) to estimate settlement-style
connectivity between the offshore wind farms in the German North Sea EEZ, for particles standing in
for e.g. bivalve larvae: released near the seabed, and settling only where they touch it.

## Quick start

OceanTracker is included as a git submodule (`oceantracker/`, pinned to a commit on its `dev`
branch), so clone with `--recurse-submodules`:

```bash
git clone --recurse-submodules <repo-url> windfarm_connectivity
cd windfarm_connectivity
# already cloned without it?  ->  git submodule update --init
```

Create the conda environment from the repo root (it installs the submodule with `pip install -e`):

```bash
mamba env create -f environment.yml             # or: conda env create -f environment.yml
mamba activate oceantracker
```

The hindcast location is set by `HINDCAST_DIR` at the top of `run_connectivity.py`
(`/hpcfreenas/hindcast/Europe/GermanBight`); change it if the hindcast lives elsewhere on your
machine. Then:

```bash
python fetch_windfarms.py                       # optional: re-download + clean polygons -> data/ (already included)
python -u run_connectivity.py --quick           # smoke test of both runs (scaled down)
python -u run_connectivity.py                   # full runs: tracks + stats (see timings below)
jupyter nbconvert --to notebook --execute --inplace windfarm_connectivity.ipynb
```

`run_connectivity.py stats` or `... tracks` runs only one of the two simulations.

`windfarm_connectivity.ipynb` contains all the analysis and plotting. It runs the simulations itself if their output is missing

## Wind farm polygons

`fetch_windfarms.py` downloads

| layer | source | licence |
|---|---|---|
| wind farm outlines | EMODnet Human Activities WFS, `emodnet:windfarmspoly`, `country='Germany'` | CC-BY 4.0 |
| German EEZ | Marine Regions (VLIZ) WFS, `MarineRegions:eez`, MRGID 5669 | CC-BY 4.0 |

and keeps North Sea farms inside the German EEZ with status *Production*, *Construction* or
*Approved* (`--include-planned` adds the large planned FEP zones/sites, many of which overlap
existing farms and get clipped). The result is **26 polygons**.

OceanTracker's polygon statistics need **single-ring, non-overlapping** polygons, so:

- *Borkum Riffgrund West*, *West II* and *OWP West* (approved) are dropped – together they are
  *Borkum Riffgrund 3* (under construction), which is kept.
- *Trianel Borkum Phase 1* and *Phase 2* are interleaved multipolygons in EMODnet and are merged
  into one *Trianel Borkum* polygon.
- Multipolygon parts < 5 % of a farm's area are dropped (e.g. a 1 km² detached piece of
  *Borkum Riffgrund 2*); remaining multipart farms (*Sandbank*) are closed with a convex hull that
  is clipped against neighbouring farms.
- Vertices are simplified to 25 m and any residual overlaps are cut away; the script asserts the
  final polygons don't overlap.

Outputs in `data/`:

- `windfarms_polygons_utm32n.json` – `{"crs", "polygons": [{"name", "points"}]}`, ready for
  `PolygonRelease` `points` and statistics `polygon_list`
- `windfarms_utm32n.geojson`, `windfarms_lonlat.geojson` – same polygons with EMODnet attributes
  (status, year, turbines, MW, area)
- `eez_germany_lonlat.geojson`, `raw/` – EEZ outline and untouched downloads

## Model set-up (`run_connectivity.py`)

Two simulations with the same physics and release polygons:

| item | value |
|---|---|
| reader | SCHISM v5 split output (auto-detected), 2020-01-01 01:00 – 2020-01-31 00:00, hourly |
| CRS | UTM 32N (EPSG:32632) – not stored in the files, inferred from the grid extent and checked against coast/EEZ and farm depths |
| releases | one `PolygonRelease` per farm, every hour from the hindcast start until 7 days before the run ends, 10 cm above the seabed (`release_at_bottom=True`, `release_offset_from_surface_or_bottom=0.1`) |
| particle lifetime | `max_age` = 14 days |
| physics | 3D velocities, `time_step` 300 s, horizontal random walk `A_H` = 1 m² s⁻¹, vertical random walk from SCHISM `diffusivity` (`use_A_Z_profile=True`) |
| **stats** run | 1000 particles per farm per pulse, no tracks; one settlement age-based polygon statistic, updated every time step |
| **tracks** run | 2 particles per farm per pulse, tracks every 6 h, used only for the maps |

**Particle numbers (stats run).** 26 farms × 1000 particles/h = 26,000 particles/h. Nothing
leaves the domain (no open boundary), and particles live up to and including 14 days (337 hourly
pulses), so the number alive grows for 14 days to 26,000 × 337 ≈ **8.76 M**. That is the maximum:
from 15 Jan settled particles are removed, and from 24 Jan no new ones are released, so fewer are
alive after that. With the current hindcast, 552 hourly pulses release about 14.4 M particles.

**Release window.** Larvae released less than 7 days before the run ends can never reach the
settlement window, so `run_connectivity.py` stops releasing then (`end` on each release group). The
time is computed at start-up from the hindcast's velocity files (`hindcast_time_range()` and
`last_useful_release()`): run end − 7 days, where the run ends at the hindcast end, or earlier with
`max_run_duration`. For the current hindcast that's 24 Jan 00:00. If more hindcast days are added,
the release end moves with them; nothing in the script needs changing. The tracks run uses the same
window, so the maps show the particles that can settle. Memory use is about 15–20 GB, mostly the
hindcast buffer.

### Settlement connectivity

The age-based statistic (`PolygonStats2D_ageBased`) uses `start=2020-01-15T01:00` (after the 14-day
spin-up), `min_age_to_bin` = 7 d, `max_age_to_bin` = 14 d, `age_bin_size` = 1 d,
`status_list=['on_bottom']`, `counting_probability=0.1`, `max_count_per_particle=1` and
`kill_when_max_counted=True`, with the release polygons as counting polygons, and `update_interval`
equal to the model time step (300 s). At every time step, a particle that is on the seabed inside a
farm while 7–14 days old settles there with a 10 % chance, so each time step on the bed is one
chance. Each particle settles **once** (possibly in its own farm) and is then removed from the
simulation.

"On the seabed" is OceanTracker's `on_bottom` status, set when a time step leaves a particle at the
bed. Within a time step the statistics are updated before resuspension lifts such particles off
again, so the statistics see them, while the screen output (written after resuspension) always shows
`Bottom:0`. In a 4-day smoke test with the earlier 900 s time step, about 9 % of the particles
inside farms were on the bed at any one update.

The notebook turns the counts into probabilities. The denominator for each source farm is the number
of its particles that are 7–14 days old at one or more updates during counting. Release
cohorts from 1–8 Jan and 17–24 Jan are only partly seen, so the probabilities are somewhat lower
than for a complete 7–14-day window (the notebook shows this in a diagram).


## Environment notes

`environment.yml` creates a conda env `oceantracker` with python 3.12, OceanTracker's dependencies
(from `oceantracker/pyproject.toml`), geopandas/shapely/pandas for the polygons and analysis, and
`ipykernel`/`nbformat`/`nbclient`/`nbconvert` for running the notebook. OceanTracker itself is
installed editable from the submodule, so the code in `oceantracker/` is what runs.

To update OceanTracker to a newer commit:

```bash
git -C oceantracker pull origin dev
git add oceantracker && git commit -m "Update oceantracker"
```

If the dependencies change, update the env with `mamba env update -f environment.yml --prune`. 