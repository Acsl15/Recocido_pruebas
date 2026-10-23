"""Inputs in SI units. Default values are demonstrations, not Ternium data."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent


def get_config(mode='demo'):
    """Return independent settings for two runs with identical boundaries."""
    if mode not in ('demo', 'tables'):
        raise ValueError("mode must be 'demo' or 'tables'")
    return {
        'mode': mode,
        'status': 'ILLUSTRATIVE - NOT TERNIUM DATA' if mode == 'demo'
                  else 'USER TABLES - INDUSTRIAL VALIDATION PENDING',
        'geometry': {'inner_diameter': 0.60, 'outer_diameter': 1.60,
                     'width': 1.00, 'n_r': 13, 'n_z': 17},
        'initial_temperature': 298.15,
        # Prescribed uniform gas: 25 -> 800 C in four hours, then hold.
        'gas_times': [0.0, 4 * 3600.0, 48 * 3600.0],
        'gas_temperatures': [298.15, 1073.15, 1073.15],
        'duration': 24 * 3600.0,
        'max_dt': 20.0,
        'safety': 0.85,
        'snapshot_interval': 1800.0,
        'targets_celsius': [580.0, 680.0, 780.0],
        'materials': [
            {'name': 'Demo A - lower diffusivity' if mode == 'demo' else 'Steel A - user table', 'rho': 7850.0,
             'cp': 600.0, 'k': 25.0, 'table': str(ROOT / 'data' / 'steel_low.csv')},
            {'name': 'Demo B - higher diffusivity' if mode == 'demo' else 'Steel B - user table', 'rho': 7850.0,
             'cp': 600.0, 'k': 45.0, 'table': str(ROOT / 'data' / 'steel_high.csv')},
        ],
        # The existing Yang implementation is the only radial model.
        # pressure/hardness below are MPa; furnace gas pressure is in Pa.
        # All values must be checked for the actual coil before industrial use.
        'radial': {'s': 0.001, 'pressure': 8.0, 'hardness': 1133.86,
                   'roughness': 3.2215e-6, 'tan_theta': 0.10,
                   'epsilon': 0.08554, 'oxide_thickness': 10e-6,
                   'k_oxide': 3.0, 'eta': 1.0},
        'boundary': {
            # Prescribed coefficients already represent the chosen total exchange.
            'mode': 'prescribed',
            'h': {'inner': 80.0, 'outer': 80.0, 'bottom': 80.0, 'top': 80.0},
            # Optional Yang (5)-(8) branch. Supply measured channel geometry,
            # actual circulating flow (not makeup flow), psi and xi to enable it.
            'q_actual': None, 'n_coils': 4, 'coil_index': 1,
            'zeta': 1 / 0.909, 'psi': None, 'xi': None,
            'channels': None,
        },
    }
