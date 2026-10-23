"""Report actual target-time sensitivity; python verify_convergence.py."""

import json
from pathlib import Path
from config import get_config
from main import run_comparison
from analysis import energy_diagnostics


def main():
    cases = [('base', 13, 17, 20.0), ('half_dt', 13, 17, 10.0),
             ('fine_mesh', 25, 33, 10.0), ('finer_mesh', 37, 49, 10.0)]
    report = {}
    for name, nr, nz, dt in cases:
        settings = get_config()
        settings['geometry'].update(n_r=nr, n_z=nz)
        settings['max_dt'] = dt
        results, rows = run_comparison(settings)
        report[name] = {'n_r': nr, 'n_z': nz, 'max_dt': dt, 'targets': rows,
                        'diagnostics': [energy_diagnostics(r) for r in results]}
    for name in ('half_dt', 'fine_mesh', 'finer_mesh'):
        reference = {'half_dt': 'base', 'fine_mesh': 'half_dt',
                     'finer_mesh': 'fine_mesh'}[name]
        report[name]['relative_changes'] = []
        for current, previous in zip(report[name]['targets'], report[reference]['targets']):
            changes = {'target_C': current['target_C'], 'reference': reference}
            for key in ('time_A_s', 'time_B_s'):
                a, b = current[key], previous[key]
                changes[key] = None if a is None or b in (None, 0) else (a-b)/b
            report[name]['relative_changes'].append(changes)
    output = Path('results_demo')
    output.mkdir(exist_ok=True)
    (output/'convergence.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v.get('relative_changes') for k, v in report.items()}, indent=2))


if __name__ == '__main__':
    main()
