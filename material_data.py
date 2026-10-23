"""Property adapters and exact enthalpy for piecewise-linear cp tables."""

from pathlib import Path
import numpy as np


def read_table(path, columns):
    """CSV without header (or a header beginning with #); SI and kelvin."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Missing {path}. Use --mode demo or supply verified tables.')
    table = np.loadtxt(path, delimiter=',', ndmin=2)
    if table.shape[0] < 2 or table.shape[1] != columns:
        raise ValueError(f'{path}: expected at least two rows and {columns} columns')
    if not np.all(np.isfinite(table)) or np.any(table <= 0):
        raise ValueError(f'{path}: all entries must be finite and positive')
    if np.any(np.diff(table[:, 0]) <= 0):
        raise ValueError(f'{path}: temperatures must be strictly increasing')
    return table


def check_range(values, lower, upper, name):
    values = np.asarray(values, dtype=float)
    if (not np.all(np.isfinite(values)) or np.any(values < lower)
            or np.any(values > upper)):
        raise ValueError(f'{name} must be finite and in [{lower}, {upper}]')
    return values


def material_from_table(table, name, source, property_function=None):
    """Build functions from T[K], rho[kg/m3], cp[J/kg/K], k[W/m/K].

    Fixed mesh and mass require constant density. Thermal expansion is excluded.
    Enthalpy is per volume (J/m3), with zero at the first tabulated temperature.
    """
    table = np.asarray(table, dtype=float)
    if (table.ndim != 2 or table.shape[1] != 4 or table.shape[0] < 2
            or not np.all(np.isfinite(table)) or np.any(table <= 0)
            or np.any(np.diff(table[:, 0]) <= 0)):
        raise ValueError('Invalid steel table')
    if not np.allclose(table[:, 1], table[0, 1], rtol=1e-10, atol=0):
        raise ValueError('Variable density requires a mass/expansion model; use constant density')
    temperatures, density, cp, k = table.T
    rho = float(density[0])
    slopes = np.diff(cp) / np.diff(temperatures)
    cumulative = np.r_[0.0, np.cumsum(0.5 * (cp[:-1] + cp[1:])
                                     * np.diff(temperatures))]

    def properties(T):
        T = check_range(T, temperatures[0], temperatures[-1], 'Steel temperature')
        if property_function is not None:
            return property_function(T)
        return (np.full_like(T, rho), np.interp(T, temperatures, cp),
                np.interp(T, temperatures, k))

    def enthalpy(T):
        T = check_range(T, temperatures[0], temperatures[-1], 'Steel temperature')
        i = np.clip(np.searchsorted(temperatures, T, side='right') - 1, 0, len(slopes)-1)
        dT = T - temperatures[i]
        return rho * (cumulative[i] + cp[i] * dT + 0.5 * slopes[i] * dT**2)

    def temperature_from_enthalpy(energy):
        e = np.asarray(energy, dtype=float) / rho
        # Only roundoff at a tabulated endpoint may be clipped.
        tolerance = 1e-10 * max(1.0, cumulative[-1])
        check_range(e, -tolerance, cumulative[-1] + tolerance, 'Specific enthalpy')
        e = np.clip(e, 0.0, cumulative[-1])
        i = np.clip(np.searchsorted(cumulative, e, side='right') - 1, 0, len(slopes)-1)
        de = e - cumulative[i]
        # Stable quadratic root; also works for zero cp slope.
        dT = 2 * de / (cp[i] + np.sqrt(cp[i]**2 + 2 * slopes[i] * de))
        return temperatures[i] + dT

    return {'name': name, 'source': source, 'properties': properties,
            'enthalpy': enthalpy, 'temperature_from_enthalpy': temperature_from_enthalpy,
            'rho_c_min': rho * float(cp.min()),
            'temperature_range': (float(temperatures[0]), float(temperatures[-1]))}


def load_models(settings):
    """Keep missing industrial tables separate from explicit demo properties."""
    if settings['mode'] == 'demo':
        materials = []
        for item in settings['materials']:
            table = [[T, item['rho'], item['cp'], item['k']] for T in (250.0, 1500.0)]
            materials.append(material_from_table(table, item['name'], 'ILLUSTRATIVE CONSTANTS'))

        def gas_properties(T):
            T = check_range(T, 250.0, 1500.0, 'Demo hydrogen temperature')
            # Deliberately simple demo values, not a CoolProp substitute.
            rho = 101325.0 * 0.002016 / (8.314462618 * T)
            return rho, np.full_like(T, 1.5e-5), np.full_like(T, 0.30), np.full_like(T, 14300.0)
        return materials, gas_properties

    # Validate both files before importing the team's eager-loading module.
    root = Path(__file__).resolve().parent
    read_table(root / 'data' / 'hydrogen.csv', 5)
    read_table(root / 'data' / 'steel.csv', 4)
    import properties
    materials = []
    for item in settings['materials']:
        path = Path(item['table']).resolve()
        table = read_table(path, 4)
        function = properties.steel_properties if path == root / 'data' / 'steel.csv' else None
        materials.append(material_from_table(table, item['name'], str(path), function))

    def gas_properties(T):
        T = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(T)):
            raise ValueError('Gas temperature must be finite')
        return properties.hydrogen_properties(T)
    return materials, gas_properties
