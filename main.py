"""Run: python main.py --mode demo --output results_demo"""

import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np

from config import get_config
from mesh import build_mesh
from material_data import load_models
from solver import solve
from analysis import compare_times, energy_diagnostics


def run_comparison(settings):
    geometry = settings['geometry']
    for key in ('n_r', 'n_z'):
        if not isinstance(geometry[key], int) or geometry[key] < 3:
            raise ValueError(f'{key} must be an integer >= 3')
    for key in ('inner_diameter', 'outer_diameter', 'width'):
        if not np.isfinite(geometry[key]) or geometry[key] <= 0:
            raise ValueError(f'Invalid {key}')
    mesh = build_mesh(r_in=geometry['inner_diameter']/2,
                      r_out=geometry['outer_diameter']/2, width=geometry['width'],
                      n_r=geometry['n_r'], n_z=geometry['n_z'], geometry='axisymmetric')
    materials, gas = load_models(settings)
    if len(materials) != 2:
        raise ValueError('Configure exactly two materials')
    results = []
    for material in materials:
        print('Simulating:', material['name'], flush=True)
        results.append(solve(mesh, material, gas, settings))
    rows = compare_times(results, settings['targets_celsius'])
    return results, rows


def save_results(settings, results, rows, output):
    from graph import save_figures
    import matplotlib
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'target_times.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    diagnostics = {}
    for i, result in enumerate(results):
        diagnostics[result['name']] = energy_diagnostics(result)
        columns = ['time', 'cold', 'hot', 'mean', 'cold_r', 'cold_z', 'hot_r', 'hot_z',
                   'energy_change', 'boundary_energy', 'residual']
        np.savetxt(output/f'history_{i+1}.csv', np.column_stack([result[c] for c in columns]),
                   delimiter=',', header=','.join(columns), comments='')
        np.savez_compressed(output/f'fields_{i+1}.npz', T=result['fields'],
                            time=result['field_times'], r=result['mesh']['r'],
                            z=result['mesh']['z'], volume_per_radian=result['mesh']['volume'],
                            status=result['status'], name=result['name'])
    metadata = {'settings': settings, 'diagnostics': diagnostics,
                'software': {'python': sys.version.split()[0], 'numpy': np.__version__,
                             'matplotlib': matplotlib.__version__},
                'target_criterion': 'first crossing of minimum nodal temperature',
                'csv_units': 'time s; temperatures K except target_C; positions m; energy J for full coil',
                'missing_crossing': 'empty CSV cell means not reached; never extrapolated'}
    (output/'run_metadata.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding='utf-8')
    save_figures(results, rows, settings['targets_celsius'], output)
    return diagnostics


def main():
    parser = argparse.ArgumentParser(description='Axisymmetric annealing demonstration; not industrial validation')
    parser.add_argument('--mode', choices=['demo', 'tables'], default='demo')
    parser.add_argument('--output', default='results_demo')
    parser.add_argument('--nr', type=int)
    parser.add_argument('--nz', type=int)
    parser.add_argument('--duration-hours', type=float)
    parser.add_argument('--max-dt', type=float)
    args = parser.parse_args()
    settings = get_config(args.mode)
    for argument, key in ((args.nr, 'n_r'), (args.nz, 'n_z')):
        if argument is not None:
            settings['geometry'][key] = argument
    if args.duration_hours is not None:
        settings['duration'] = args.duration_hours*3600
    if args.max_dt is not None:
        settings['max_dt'] = args.max_dt
    print(settings['status'], flush=True)
    results, rows = run_comparison(settings)
    diagnostics = save_results(settings, results, rows, args.output)
    for row in rows:
        print(row)
    for name, values in diagnostics.items():
        print(name, 'energy residual:', f"{values['relative_residual']:.3e}")
    print('Saved:', Path(args.output).resolve())


if __name__ == '__main__':
    main()
