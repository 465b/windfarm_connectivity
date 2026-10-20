"""
OceanTracker run: connectivity between offshore wind farms in the German EEZ, using the
SCHISM (v5 split output) GermanBight hindcast (2020-01-01 01:00 to 2020-01-31 00:00).

Particles stand in for e.g. bivalve larvae: released 10 cm above the seabed, passive (no swimming
or sinking), and they can only settle while touching the seabed.

Two runs, same physics and release polygons:
  stats   one PolygonRelease per wind farm, 1000 particles every hour until 7 days before the run
          ends (computed from the hindcast files, so it follows an extended hindcast): later larvae
          can't reach 7 days, the start of the settlement window. max_age 14 days
          (at most ~8.76 M particles alive, see below).
          On-the-fly statistics only, no tracks:
            settlement-style age-based polygon stats: from day 14, a particle on the seabed
            (status on_bottom) inside a farm while aged 7-14 days settles with a 10 % probability
            at each update, which is every model time step (300 s); a settled particle is removed
            from the simulation (killed), so it can't be counted again
            (PolygonStats2D_ageBased with status_list=['on_bottom'], counting_probability=0.1,
            max_count_per_particle=1, kill_when_max_counted=True, update_interval=time_step).
  tracks  same releases with 2 particles per pulse, tracks written every 6 h, for plotting

Particles in the system: 26 farms x 1000 particles/h = 26,000 particles/h. Nothing leaves the
domain (no open boundary), and particles live up to and including max_age = 14 d, i.e. 337 hourly
pulses, so the number alive levels off at 26,000 x 337 = 8.76 M on 2020-01-15. That is the maximum:
from 2020-01-15 settled particles are removed, and from the release end (2020-01-24 00:00 for the
current 30-day hindcast) no new ones are added, so fewer are alive after that. With the current
hindcast, 552 hourly pulses give about 14.4 M particles released in total.

Usage (use python -u when redirecting to a log file, otherwise progress output is block-buffered)
  python run_connectivity.py                 # stats + tracks runs
  python run_connectivity.py stats           # only the statistics run
  python run_connectivity.py tracks          # only the tracks run
  python run_connectivity.py --quick stats   # short smoke test (scaled-down ages and particles)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr
from oceantracker.main import OceanTracker

HINDCAST_DIR = '/hpcfreenas/hindcast/Europe/GermanBight'
VELOCITY_FILE_MASK = 'horizontalVelX_*.nc'   # OceanTracker takes the hindcast time axis from the velocity files
HERE = Path(__file__).parent
TIME_STEP = 300.   # model time step [s], also the settlement statistic's update interval

hour = 3600.
day = 24 * hour


def hindcast_time_range():
    '''First and last time step of the hindcast, as numpy datetime64.'''
    times = []
    for f in Path(HINDCAST_DIR).glob(VELOCITY_FILE_MASK):
        with xr.open_dataset(f) as ds:
            times.append(ds['time'].values[[0, -1]])
    times = np.concatenate(times)
    return times.min(), times.max()


def last_useful_release(min_age_to_count, max_run_duration=None):
    '''Time of the last release whose particles can still reach min_age_to_count before the run ends.
    Releases start at the hindcast start, and the run ends at the hindcast end or after max_run_duration,
    so this follows the hindcast when it is extended.'''
    start, end = hindcast_time_range()
    if max_run_duration is not None:
        end = min(end, start + np.timedelta64(int(max_run_duration), 's'))
    return str((end - np.timedelta64(int(min_age_to_count), 's')).astype('datetime64[s]'))


def _base_model(run_name, write_tracks, pulse_size, max_age, release_end, max_run_duration=None, processors=None):
    with open(HERE / 'data' / 'windfarms_polygons_utm32n.json') as f:
        polygons = json.load(f)['polygons']

    ot = OceanTracker()
    settings = dict(
        run_output_dir=str(HERE / 'output' / run_name),
        time_step=TIME_STEP,
        use_A_Z_profile=True,        # vertical random walk from SCHISM diffusivity
        use_open_boundary=False,     # no hgrid.gr3 available, open boundary nodes are unknown
        write_tracks=write_tracks,
        screen_output_time_interval=6 * hour,
        processors=processors,
    )
    if max_run_duration is not None:
        settings['max_run_duration'] = max_run_duration
    ot.settings(**settings)
    ot.add_class('reader',
                 input_dir=HINDCAST_DIR,
                 file_mask='*.nc',
                 EPSG_code=32632)    # not used unless use_geographic_coords, recorded for reference

    ot.add_class('dispersion', A_H=1.0)  # horizontal eddy diffusivity m^2/s, plausible for km-scale mesh
                                         # vertical is dynamic based on schism output

    # release every hour from the start of the hindcast until release_end
    for p in polygons:
        ot.add_class('release_groups',
                     name=p['name'],
                     class_name='PolygonRelease',
                     points=p['points'],
                     pulse_size=pulse_size,
                     release_interval=hour,
                     end=release_end,
                     max_age=max_age,
                     release_at_bottom=True,                     # e.g. bivalve larvae near the seabed ...
                     release_offset_from_surface_or_bottom=0.1,  # ... 10 cm above it
                     )
    return ot


def build_stats_model(quick=False, pulse_size=1000, processors=None):
    if quick:
        # scaled down: ages 1-2 d counted from day 2, 4 day run
        max_age, min_age_bin, age_bin, count_start, duration = 2 * day, 1 * day, 6 * hour, '2020-01-03T01:00:00', 4 * day
        pulse_size = min(pulse_size, 20)
    else:
        max_age, min_age_bin, age_bin, count_start, duration = 14 * day, 7 * day, 1 * day, '2020-01-15T01:00:00', None

    # stop releasing once new larvae can no longer reach the start of the counting window (min_age_bin)
    # before the run ends; they could never settle
    release_end = last_useful_release(min_age_bin, duration)

    ot = _base_model('windfarm_connectivity_quick' if quick else 'windfarm_connectivity', write_tracks=False,
                     pulse_size=pulse_size, max_age=max_age, release_end=release_end,
                     max_run_duration=duration, processors=processors)

    # settlement-style farm-to-farm connectivity, by age at settlement; use_release_group_polygons
    # makes the counting polygons the same as the release polygons (same order).
    # Only particles touching the seabed can be counted ("on_bottom" is set when a particle reaches
    # the bed, and is visible to the statistics before resuspension lifts it off again). The update
    # interval is the model time step, so every time step a particle spends on the seabed gives it
    # one 10 % chance to settle. Counting starts once the age structure is spun up (after max_age).
    ot.add_class('particle_statistics',
                 name='windfarm_connectivity',
                 class_name='PolygonStats2D_ageBased',
                 use_release_group_polygons=True,
                 update_interval=TIME_STEP,    # check every time step
                 start=count_start,
                 min_age_to_bin=min_age_bin,
                 max_age_to_bin=max_age,
                 age_bin_size=age_bin,
                 status_list=['on_bottom'],    # only particles on the seabed ...
                 counting_probability=0.1,     # ... settle with a 10 % chance per time step on the bed
                 max_count_per_particle=1,     # count each particle once ...
                 kill_when_max_counted=True,   # ... then it has settled: remove it from the simulation,
                                               # so it is not counted by any other statistic
                 )

    return ot


def build_tracks_model(quick=False, pulse_size=2, processors=None):
    # same release window as the stats run, so the maps show the particles that can settle
    min_age_bin, duration = (1 * day, 4 * day) if quick else (7 * day, None)
    ot = _base_model('windfarm_tracks_quick' if quick else 'windfarm_tracks', write_tracks=True,
                     pulse_size=pulse_size, max_age=2 * day if quick else 14 * day,
                     release_end=last_useful_release(min_age_bin, duration),
                     max_run_duration=duration, processors=processors)
    ot.add_class('tracks_writer', update_interval=6 * hour)
    return ot


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('runs', nargs='?', default='both', choices=['both', 'stats', 'tracks'])
    parser.add_argument('--quick', action='store_true', help='short smoke test')
    parser.add_argument('--pulse-size', type=int, default=1000, help='stats run: particles per farm per hourly release')
    parser.add_argument('--processors', type=int, default=None, help='max threads (default: all physical cores)')
    args = parser.parse_args()

    if args.runs in ('both', 'tracks'):
        print('case info file:', build_tracks_model(quick=args.quick, processors=args.processors).run())
    if args.runs in ('both', 'stats'):
        print('case info file:', build_stats_model(quick=args.quick, pulse_size=args.pulse_size,
                                                   processors=args.processors).run())
