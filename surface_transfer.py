"""Surface coefficients; optional Yang (5)-(8), without editing team files."""

import numpy as np
from flow import surface_flows
from convection import convection_heat_flux

SURFACES = ('inner', 'outer', 'bottom', 'top')


def positive(value, name, allow_zero=False):
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'{name} requires an explicit numerical value') from error
    if not np.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        raise ValueError(f'{name} must be finite and positive (zero allowed: {allow_zero})')
    return value


def surface_coefficients(T_g, boundary, gas_properties):
    """Return h_eff in W/m2/K for all four surfaces.

    q_actual is circulating m3/s at the gas operating conditions, not Nm3/s.
    The Yang branch uses supplied channel area, hydraulic diameter and length.
    psi and xi are explicit empirical inputs; no separate surface radiation is added.
    """
    positive(T_g, 'T_g')
    if boundary['mode'] == 'prescribed':
        return {s: positive(boundary['h'][s], f'h_{s}', True) for s in SURFACES}
    if boundary['mode'] != 'yang':
        raise ValueError("boundary mode must be 'prescribed' or 'yang'")
    n = boundary['n_coils']
    index = boundary['coil_index']
    if not isinstance(n, int) or n < 2 or not isinstance(index, int) or not 0 <= index < n:
        raise ValueError('Invalid stack size or coil_index')
    psi = positive(boundary['psi'], 'psi')
    xi = positive(boundary['xi'], 'xi')
    zeta = positive(boundary['zeta'], 'zeta')
    q = positive(boundary['q_actual'], 'q_actual')
    if not isinstance(boundary.get('channels'), dict):
        raise ValueError('Yang mode requires channel geometry for each exposed surface')
    flows = surface_flows(q_total=q, n_coils=n, zeta=zeta)[index]
    rho, mu, k, _ = gas_properties(T_g)
    rho, mu, k = [positive(v, label) for v, label in zip((rho, mu, k), ('rho_g', 'mu_g', 'k_g'))]
    result = {}
    for surface in SURFACES:
        if flows[surface] == 0:
            # Explicit model choice: no convective channel on this face.
            result[surface] = 0.0
            continue
        channel = boundary['channels'][surface]
        area = positive(channel['area'], 'flow area')
        D = positive(channel['diameter'], 'hydraulic diameter')
        L = positive(channel['length'], 'flow length')
        velocity = flows[surface] / area
        reynolds = rho * velocity * D / mu
        if reynolds < 10000:
            raise ValueError(f'{surface}: Re={reynolds:.0f}; turbulent Yang correlation not justified')
        h = 0.023 * psi * velocity**0.8 / D**0.2 * k * (mu/rho)**-0.8 * (1+(D/L)**0.7)
        result[surface] = h * (xi if T_g >= 800.0 else 1.0)
    return result


def boundary_power(mesh, T, T_g, h):
    """Nodal W per radian; call the existing scalar convection function."""
    power = np.zeros_like(T)
    for surface in SURFACES:
        area = mesh['area_' + surface]
        active = area > 0
        flux = np.array([convection_heat_flux(h_eff=h[surface], T_g=float(T_g),
                                             T_surface=float(value)) for value in T[active]])
        power[active] += flux * area[active]
    return power
