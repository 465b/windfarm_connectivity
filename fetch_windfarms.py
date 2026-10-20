"""
Fetch offshore wind farm outlines in the German EEZ (North Sea) and turn them
into non-overlapping single-ring polygons usable by OceanTracker for both
PolygonRelease release groups and polygon statistics.

Sources
  - Wind farms: EMODnet Human Activities WFS, layer emodnet:windfarmspoly
  - German EEZ: Marine Regions (VLIZ) WFS, layer MarineRegions:eez, MRGID 5669

Outputs (in ./data)
  - windfarms_lonlat.geojson      cleaned polygons, WGS84, with attributes
  - windfarms_utm32n.geojson      same, in the hindcast CRS (EPSG:32632)
  - windfarms_polygons_utm32n.json  list of {"name", "points"} for OceanTracker
  - eez_germany_lonlat.geojson    EEZ outline for plotting
  - raw/                          untouched WFS downloads

Usage
  python fetch_windfarms.py [--include-planned]
"""
import argparse
import json
import itertools
from pathlib import Path

import requests
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon

HINDCAST_EPSG = 32632  # SCHISM GermanBight hindcast is in UTM 32N (no CRS in files, checked by extent)

EMODNET_WFS = 'https://ows.emodnet-humanactivities.eu/wfs'
MARINEREGIONS_WFS = 'https://geo.vliz.be/geoserver/MarineRegions/wfs'

# EMODnet lists the three approved sub-sites that together form "Borkum Riffgrund 3"
# (under construction) as separate, fully overlapping polygons. Keep BR3 only.
DROP_NAMES = ['Borkum Riffgrund West', 'Borkum Riffgrund West II', 'OWP West']

# The two Trianel Borkum phases are interleaved multipolygons in EMODnet and cannot be
# separated into non-overlapping single rings, so treat the site as one farm.
MERGE = {'Trianel Windpark Borkum': ['Trianel Windpark Borkum Phase 1', 'Trianel Windpark Borkum Phase 2']}

MIN_PART_FRACTION = 0.05  # drop multipolygon parts smaller than this fraction of the farm area
SIMPLIFY_TOL_M = 25.      # vertex simplification tolerance, small vs. hindcast grid spacing
NORTH_SEA_MAX_LON = 9.5   # German EEZ also covers the Baltic, which is outside the hindcast


def wfs_get(url, type_name, cql_filter, out_file):
    params = dict(service='WFS', version='2.0.0', request='GetFeature', typeNames=type_name,
                  outputFormat='application/json', srsName='EPSG:4326', CQL_FILTER=cql_filter)
    r = requests.get(url, params=params, timeout=180)
    r.raise_for_status()
    out_file.write_bytes(r.content)
    return gpd.read_file(out_file)


def to_single_ring(geom, others):
    '''Reduce a (multi)polygon to one exterior ring that does not overlap "others".'''
    parts = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    total = sum(p.area for p in parts)
    parts = [p for p in parts if p.area >= MIN_PART_FRACTION * total]
    if len(parts) == 1:
        return Polygon(parts[0].exterior)
    # several substantial parts separated by narrow gaps: close the gap with the hull,
    # then cut away anything belonging to neighbouring farms
    hull = MultiPolygon(parts).convex_hull.difference(others)
    if isinstance(hull, MultiPolygon):
        hull = max(hull.geoms, key=lambda p: p.area)
    return Polygon(hull.exterior)


def main(include_planned=False):
    data_dir = Path(__file__).parent / 'data'
    raw_dir = data_dir / 'raw'
    raw_dir.mkdir(parents=True, exist_ok=True)

    wf = wfs_get(EMODNET_WFS, 'emodnet:windfarmspoly', "country='Germany'",
                 raw_dir / 'emodnet_windfarmspoly_germany.geojson')
    eez = wfs_get(MARINEREGIONS_WFS, 'MarineRegions:eez', "mrgid=5669",
                  raw_dir / 'marineregions_eez_germany.geojson')
    print(f'downloaded {len(wf)} German wind farm polygons')

    status = ['Production', 'Construction', 'Approved'] + (['Planned'] if include_planned else [])
    wf = wf[wf['status'].isin(status) & ~wf['name'].isin(DROP_NAMES)]
    wf = wf[wf.geometry.representative_point().x < NORTH_SEA_MAX_LON]
    wf = wf[wf.geometry.representative_point().within(eez.geometry.iloc[0])]
    wf = wf.to_crs(HINDCAST_EPSG)

    for new_name, names in MERGE.items():
        sel = wf['name'].isin(names)
        if sel.any():
            merged = wf[sel].iloc[[0]].copy()
            merged['name'] = new_name
            merged[['n_turbines', 'power_mw']] = wf.loc[sel, ['n_turbines', 'power_mw']].sum().values
            # +/-100 m buffer closes the slits left between the interleaved phases
            merged = merged.set_geometry([wf[sel].geometry.union_all().buffer(100.).buffer(-100.)],
                                         crs=HINDCAST_EPSG)
            wf = gpd.GeoDataFrame(pd.concat([wf[~sel], merged]), crs=HINDCAST_EPSG)
    wf = wf.sort_values('name').reset_index(drop=True)

    # single exterior ring per farm, no overlaps
    # (earlier rows win where polygons touch/overlap; for operating farms overlaps are ~0)
    geoms = []
    for n, g in enumerate(wf.geometry):
        others = gpd.GeoSeries(geoms + list(wf.geometry[n + 1:])).union_all()
        p = to_single_ring(g, others).simplify(SIMPLIFY_TOL_M, preserve_topology=True)
        if geoms:
            p = p.difference(gpd.GeoSeries(geoms).union_all())
            if isinstance(p, MultiPolygon):
                p = max(p.geoms, key=lambda q: q.area)
        geoms.append(Polygon(p.exterior))
    wf = wf.set_geometry(gpd.GeoSeries(geoms, crs=HINDCAST_EPSG))

    for (i, a), (j, b) in itertools.combinations(enumerate(wf.geometry), 2):
        assert a.intersection(b).area < 1., f'polygons overlap: {wf.name[i]} / {wf.name[j]}'
    assert all(g.is_valid and isinstance(g, Polygon) for g in wf.geometry)

    keep = ['name', 'status', 'year', 'n_turbines', 'power_mw', 'type_inst', 'dist_coast']
    wf = wf[keep + ['geometry']]
    wf['area_km2'] = wf.geometry.area / 1e6
    wf['short_name'] = wf['name'].str.replace(r'\s*\(.*\)', '', regex=True).str.replace('Trianel Windpark ', 'Trianel ')

    wf.to_file(data_dir / 'windfarms_utm32n.geojson', driver='GeoJSON')
    wf.to_crs(4326).to_file(data_dir / 'windfarms_lonlat.geojson', driver='GeoJSON')
    eez.to_file(data_dir / 'eez_germany_lonlat.geojson', driver='GeoJSON')

    polygon_list = [dict(name=r.short_name, points=[list(xy) for xy in r.geometry.exterior.coords])
                    for r in wf.itertuples()]
    with open(data_dir / 'windfarms_polygons_utm32n.json', 'w') as f:
        json.dump(dict(crs=f'EPSG:{HINDCAST_EPSG}', polygons=polygon_list), f, indent=1)

    print(wf[['short_name', 'status', 'power_mw', 'area_km2']].to_string())
    print(f'wrote {len(wf)} polygons to {data_dir}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--include-planned', action='store_true',
                        help='also include planned sites/zones (many overlap existing farms and are clipped)')
    main(parser.parse_args().include_planned)
