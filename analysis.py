"""Cold/hot histories, first crossings and comparison of two materials."""

import numpy as np


def target_time(times, temperatures, target):
    """First upward crossing, linearly interpolated; None if never reached.

    If initially at/above the target return the initial time. No extrapolation.
    Use the cold-point history to require every mesh node to reach the target.
    """
    times, temperatures = np.asarray(times, dtype=float), np.asarray(temperatures, dtype=float)
    if (times.ndim != 1 or times.size == 0 or temperatures.shape != times.shape
            or not np.all(np.isfinite(times)) or not np.all(np.isfinite(temperatures))
            or not np.isfinite(target) or np.any(np.diff(times) <= 0)):
        raise ValueError('Invalid time/temperature history or target')
    reached = np.flatnonzero(temperatures >= target)
    if reached.size == 0:
        return None
    i = int(reached[0])
    if i == 0:
        return float(times[0])
    fraction = (target-temperatures[i-1])/(temperatures[i]-temperatures[i-1])
    return float(times[i-1] + fraction*(times[i]-times[i-1]))


def compare_times(results, targets_celsius):
    """Difference is A minus B; ratios undefined when a target is not reached."""
    if len(results) != 2:
        raise ValueError('Exactly two materials are required')
    rows = []
    for target in targets_celsius:
        if not np.isfinite(target) or not 580 <= target <= 780:
            raise ValueError('Project targets must be between 580 and 780 C')
        a, b = [target_time(r['time'], r['cold'], target+273.15) for r in results]
        rows.append({'target_C': float(target), 'time_A_s': a, 'time_B_s': b,
                     'difference_A_minus_B_s': None if a is None or b is None else a-b})
    return rows


def energy_diagnostics(result):
    """Whole-coil joules; normalize by maximum cumulative boundary energy."""
    error = float(np.max(np.abs(result['residual'])))
    scale = max(1.0, float(np.max(np.abs(result['boundary_energy']))))
    return {'maximum_absolute_residual_J': error, 'relative_residual': error/scale,
            'steps': result['steps'], 'minimum_dt_s': result['dt_min'],
            'maximum_dt_s': result['dt_max'],
            'maximum_dt_over_safe_limit': result['max_step_ratio'],
            'alpha_r_range_m2_s': result['alpha_r_range'],
            'alpha_z_range_m2_s': result['alpha_z_range']}
