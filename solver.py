"""Conservative explicit heat balance on the team's axisymmetric mesh."""

import numpy as np
from mesh import stable_time_step
from radial_conductivity import effective_radial_conductivity
from surface_transfer import surface_coefficients, boundary_power, positive


def face_conductances(mesh, k_r, k_z):
    """Harmonic interface conductivity, using mesh areas and spacings."""
    shape = mesh['volume'].shape
    k_r = np.broadcast_to(np.asarray(k_r, dtype=float), shape)
    k_z = np.broadcast_to(np.asarray(k_z, dtype=float), shape)
    if any(np.any(~np.isfinite(k)) or np.any(k <= 0) for k in (k_r, k_z)):
        raise ValueError('Conductivities must be finite and positive')
    radial = 2 / (1/k_r[:, :-1] + 1/k_r[:, 1:])
    axial = 2 / (1/k_z[:-1, :] + 1/k_z[1:, :])
    return radial * mesh['area_r'] / mesh['dr'], axial * mesh['area_z'] / mesh['dz']


def conduction_power(T, G_r, G_z):
    """Each internal face adds equal and opposite powers to its two cells."""
    power = np.zeros_like(T)
    radial = G_r * (T[:, 1:] - T[:, :-1])
    axial = G_z * (T[1:, :] - T[:-1, :])
    power[:, :-1] += radial
    power[:, 1:] -= radial
    power[:-1, :] += axial
    power[1:, :] -= axial
    return power


def event_times(settings):
    """Land exactly on gas knots, output times and Yang's 800 K transition."""
    times = np.asarray(settings['gas_times'], dtype=float)
    temperatures = np.asarray(settings['gas_temperatures'], dtype=float)
    end = settings['duration']
    snapshots = np.r_[np.arange(0, end, settings['snapshot_interval']), end]
    events = list(times) + list(snapshots)
    for a, b, Ta, Tb in zip(times[:-1], times[1:], temperatures[:-1], temperatures[1:]):
        if (Ta < 800 < Tb) or (Tb < 800 < Ta):
            events.append(a + (800-Ta) / (Tb-Ta) * (b-a))
    events = np.unique([t for t in events if 0 < t <= end])
    return events, snapshots


def validate_settings(mesh, settings):
    if mesh['geometry'] != 'axisymmetric':
        raise ValueError('This solver requires the axisymmetric mesh (per radian)')
    for key in ('duration', 'max_dt', 'snapshot_interval'):
        positive(settings[key], key)
    if not 0 < settings['safety'] <= 1:
        raise ValueError('safety must be in (0, 1]')
    times = np.asarray(settings['gas_times'], dtype=float)
    temperatures = np.asarray(settings['gas_temperatures'], dtype=float)
    if (times.ndim != 1 or times.size < 2 or temperatures.shape != times.shape
            or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0)
            or times[0] != 0 or times[-1] < settings['duration']
            or not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0)):
        raise ValueError('Invalid gas profile; it must cover [0, duration] in seconds and kelvin')


def solve(mesh, material, gas_properties, settings):
    """Solve one material, retaining every-step diagnostics and sparse 2D maps.

    Advance volumetric enthalpy, not cp(T)*T: temperature-dependent cp then
    conserves energy. The global minimum rho*cp makes the explicit step bound
    safe even when a step crosses a decreasing-cp table segment.
    No internal sources, expansion or gas energy equation are included.
    """
    validate_settings(mesh, settings)
    volume = mesh['volume']
    T = np.broadcast_to(np.asarray(settings['initial_temperature'], dtype=float), volume.shape).copy()
    U = material['enthalpy'](T)
    initial_energy = float(np.sum(U * volume))
    rho_c_min = positive(material['rho_c_min'], 'minimum rho*cp')
    events, snapshots = event_times(settings)
    next_event = 0
    next_snapshot = 1
    time = 0.0
    boundary_energy = 0.0
    steps = 0
    dt_min, dt_max = np.inf, 0.0
    max_step_ratio = 0.0
    alpha_r_min, alpha_r_max = np.inf, 0.0
    alpha_z_min, alpha_z_max = np.inf, 0.0
    history = {key: [] for key in ('time', 'cold', 'hot', 'mean', 'cold_r', 'cold_z',
                                  'hot_r', 'hot_z', 'energy_change', 'boundary_energy', 'residual')}
    fields = [T.copy()]
    field_times = [0.0]

    def record():
        cold_index = np.unravel_index(np.argmin(T), T.shape)
        hot_index = np.unravel_index(np.argmax(T), T.shape)
        energy_change = float(np.sum(U * volume)) - initial_energy
        values = (time, T[cold_index], T[hot_index], np.sum(T*volume)/np.sum(volume),
                  mesh['r'][cold_index[1]], mesh['z'][cold_index[0]],
                  mesh['r'][hot_index[1]], mesh['z'][hot_index[0]],
                  energy_change*2*np.pi, boundary_energy*2*np.pi,
                  (energy_change-boundary_energy)*2*np.pi)
        for key, value in zip(history, values):
            history[key].append(float(value))

    record()
    while time < settings['duration'] - 1e-9:
        T_g = float(np.interp(time, settings['gas_times'], settings['gas_temperatures']))
        rho, cp, k_z = [np.broadcast_to(np.asarray(v, dtype=float), T.shape)
                        for v in material['properties'](T)]
        if any(np.any(~np.isfinite(v)) or np.any(v <= 0) for v in (rho, cp, k_z)):
            raise ValueError('Invalid steel properties')
        if np.any(rho*cp < rho_c_min * (1-1e-12)):
            raise ValueError('rho_c_min is not a lower bound for this material')
        # The interlayer gas conductivity is evaluated at local steel temperature;
        # the uniform furnace temperature selects Yang's low/high branch.
        k_g = gas_properties(T)[2]
        k_r = effective_radial_conductivity(T_m=T, T_g=T_g, k_s=k_z, k_g=k_g,
                                            **settings['radial'])
        h = surface_coefficients(T_g, settings['boundary'], gas_properties)
        G_r, G_z = face_conductances(mesh, k_r, k_z)
        limit = stable_time_step(mesh=mesh, rho_c=np.full_like(T, rho_c_min),
                                 G_r=G_r, G_z=G_z, h=h, safety=settings['safety'])
        dt = min(settings['max_dt'], limit, events[next_event]-time)
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Nonpositive or invalid time step')
        outside = boundary_power(mesh, T, T_g, h)
        power = conduction_power(T, G_r, G_z) + outside
        U_new = U + dt * power / volume
        T_new = material['temperature_from_enthalpy'](U_new)
        # A positive explicit balance cannot exceed the old field and reservoir.
        if (np.min(T_new) < min(float(T.min()), T_g)-1e-7
                or np.max(T_new) > max(float(T.max()), T_g)+1e-7):
            raise RuntimeError('Discrete temperature bounds violated')
        boundary_energy += dt * float(outside.sum())
        time += dt
        if abs(time-events[next_event]) < 1e-8:
            time = float(events[next_event])
            next_event += 1
        T, U = T_new, U_new
        steps += 1
        dt_min, dt_max = min(dt_min, dt), max(dt_max, dt)
        max_step_ratio = max(max_step_ratio, dt / limit)
        alpha_r_min = min(alpha_r_min, float(np.min(k_r/(rho*cp))))
        alpha_r_max = max(alpha_r_max, float(np.max(k_r/(rho*cp))))
        alpha_z_min = min(alpha_z_min, float(np.min(k_z/(rho*cp))))
        alpha_z_max = max(alpha_z_max, float(np.max(k_z/(rho*cp))))
        record()
        if next_snapshot < len(snapshots) and abs(time-snapshots[next_snapshot]) < 1e-7:
            fields.append(T.copy())
            field_times.append(time)
            next_snapshot += 1
        if steps >= 2_000_000:
            raise RuntimeError('Too many steps; check mesh, conductivities and units')

    result = {key: np.array(value) for key, value in history.items()}
    result.update({'name': material['name'], 'source': material['source'],
                   'status': settings['status'], 'fields': np.array(fields),
                   'field_times': np.array(field_times), 'mesh': mesh,
                   'steps': steps, 'dt_min': dt_min, 'dt_max': dt_max,
                   'max_step_ratio': max_step_ratio,
                   'alpha_r_range': [alpha_r_min, alpha_r_max],
                   'alpha_z_range': [alpha_z_min, alpha_z_max]})
    return result
