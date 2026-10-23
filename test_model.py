"""Independent checks: python -m unittest -v test_model"""

import copy
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

from config import get_config
from mesh import build_mesh, stable_time_step
from flow import flow_distribution_coefficients, channel_flows
from material_data import load_models, material_from_table, read_table
from radial_conductivity import effective_radial_conductivity
from surface_transfer import surface_coefficients, boundary_power
from solver import solve, face_conductances, conduction_power
from analysis import target_time, compare_times, energy_diagnostics


def small_case():
    settings = get_config()
    settings.update(duration=600.0, max_dt=10.0, snapshot_interval=100.0,
                    gas_times=[0.0, 600.0], gas_temperatures=[900.0, 900.0])
    grid = build_mesh(r_in=0.3, r_out=0.8, width=1, n_r=5, n_z=7)
    materials, gas = load_models(settings)
    return settings, grid, materials[0], gas


class PhysicalChecks(unittest.TestCase):
    def test_exact_cylinder_volume_and_surface_areas(self):
        _, grid, _, _ = small_case()
        self.assertAlmostEqual(2*np.pi*grid['volume'].sum(), np.pi*(0.8**2-0.3**2))
        self.assertAlmostEqual(2*np.pi*grid['area_inner'].sum(), 2*np.pi*0.3)
        self.assertAlmostEqual(2*np.pi*grid['area_outer'].sum(), 2*np.pi*0.8)
        for side in ('bottom', 'top'):
            self.assertAlmostEqual(2*np.pi*grid['area_'+side].sum(), np.pi*(0.8**2-0.3**2))

    def test_yang_flow_table_and_continuity(self):
        fd = flow_distribution_coefficients(n_coils=4, zeta=1/0.909)
        plates, sides, top = channel_flows(q_total=1, fd=fd)
        np.testing.assert_allclose(100*plates, [20.2, 14.6, 10.9], atol=0.05)
        np.testing.assert_allclose(100*sides, [100, 79.8, 65.2, 54.3], atol=0.05)
        self.assertAlmostEqual(plates.sum()+top, 1)

    def test_internal_exchange_conserves_energy(self):
        _, grid, _, _ = small_case()
        rng = np.random.default_rng(7)
        T = rng.uniform(300, 900, grid['volume'].shape)
        G_r, G_z = face_conductances(grid, rng.uniform(1, 20, T.shape), 40)
        self.assertLess(abs(conduction_power(T, G_r, G_z).sum()), 1e-9)

    def test_surface_power_matches_whole_cylinder(self):
        _, grid, _, _ = small_case()
        T = np.full_like(grid['volume'], 300)
        h = dict.fromkeys(('inner', 'outer', 'bottom', 'top'), 100.0)
        power = boundary_power(grid, T, 800.0, h).sum()*2*np.pi
        area = 2*np.pi*(0.3+0.8) + 2*np.pi*(0.8**2-0.3**2)
        self.assertAlmostEqual(power, 100*500*area)

    def test_enthalpy_against_analytic_linear_cp(self):
        # cp = 100 + T, exact integral is 100*(T-250)+(T^2-250^2)/2.
        material = material_from_table([[250, 7800, 350, 20], [1500, 7800, 1600, 20]], 'test', 'analytic')
        T = np.array([[250.0, 300.0], [900.0, 1500.0]])
        exact = 7800*(100*(T-250)+0.5*(T**2-250**2))
        np.testing.assert_allclose(material['enthalpy'](T), exact, rtol=1e-13)
        np.testing.assert_allclose(material['temperature_from_enthalpy'](exact), T, atol=1e-10)

    def test_enthalpy_inverse_with_falling_cp(self):
        table = [[250, 7800, 800, 25], [700, 7800, 1200, 25], [1100, 7800, 500, 25]]
        material = material_from_table(table, 'test', 'test')
        T = np.linspace(250, 1100, 101)
        np.testing.assert_allclose(material['temperature_from_enthalpy'](material['enthalpy'](T)), T, atol=1e-10)

    def test_equilibrium_remains_uniform(self):
        settings, grid, material, gas = small_case()
        settings['initial_temperature'] = 900.0
        result = solve(grid, material, gas, settings)
        np.testing.assert_allclose(result['fields'], 900, atol=1e-10)
        self.assertLess(np.max(np.abs(result['boundary_energy'])), 1e-6)

    def test_heating_bounds_and_energy(self):
        settings, grid, material, gas = small_case()
        settings['max_dt'] = 1e9  # Stability bound must override this.
        result = solve(grid, material, gas, settings)
        self.assertGreater(result['cold'][-1], settings['initial_temperature'])
        self.assertLessEqual(result['hot'].max(), 900.0)
        self.assertLessEqual(result['max_step_ratio'], 1+1e-12)
        final = result['fields'][-1]
        # Independent constant-cp energy, not the solver's energy accumulator.
        exact = 7850*600*np.sum((final-settings['initial_temperature'])*grid['volume'])*2*np.pi
        self.assertAlmostEqual(exact/result['boundary_energy'][-1], 1.0, places=10)
        self.assertLess(energy_diagnostics(result)['relative_residual'], 1e-10)

    def test_cooling_bounds_and_negative_boundary_energy(self):
        settings, grid, material, gas = small_case()
        settings['initial_temperature'] = 1000.0
        settings['gas_temperatures'] = [300.0, 300.0]
        result = solve(grid, material, gas, settings)
        self.assertGreaterEqual(result['cold'].min(), 300)
        self.assertLess(result['hot'][-1], 1000)
        self.assertLess(result['boundary_energy'][-1], 0)

    def test_insulated_variable_cp_energy(self):
        settings, grid, _, gas = small_case()
        settings['boundary']['h'] = dict.fromkeys(('inner', 'outer', 'bottom', 'top'), 0.0)
        settings['initial_temperature'] = 500+200*grid['Z']
        material = material_from_table([[250, 7800, 350, 20], [1500, 7800, 1600, 20]], 'test', 'analytic')
        result = solve(grid, material, gas, settings)
        initial, final = settings['initial_temperature'], result['fields'][-1]
        independent_change = 7800*np.sum((100*(final-initial)+0.5*(final**2-initial**2))*grid['volume'])
        self.assertLess(abs(independent_change), 1e-5)
        self.assertEqual(result['boundary_energy'][-1], 0.0)

    def test_analytic_axial_cosine_and_mesh_refinement(self):
        errors = []
        for nz in (9, 17):
            settings, _, material, gas = small_case()
            grid = build_mesh(r_in=0.3, r_out=0.8, width=1, n_r=5, n_z=nz)
            settings.update(duration=3000.0, gas_times=[0, 3000], max_dt=1.0,
                            snapshot_interval=3000.0)
            settings['boundary']['h'] = dict.fromkeys(('inner', 'outer', 'bottom', 'top'), 0.0)
            settings['initial_temperature'] = 600 + 100*np.cos(np.pi*grid['Z'])
            result = solve(grid, material, gas, settings)
            exact = 600+100*np.cos(np.pi*grid['Z'])*np.exp(-(25/(7850*600))*np.pi**2*3000)
            errors.append(float(np.max(np.abs(result['fields'][-1]-exact))))
        self.assertLess(errors[1], 0.06)
        self.assertLess(errors[1], 0.3*errors[0])

    def test_radial_uses_existing_yang_function(self):
        settings, _, _, _ = small_case()
        parameters = settings['radial'].copy()
        # Limit of no contact, no emission, no oxide: two series resistances.
        parameters.update(pressure=0, epsilon=0, oxide_thickness=0)
        k = effective_radial_conductivity(T_m=900, T_g=900, k_s=30, k_g=0.3, **parameters)
        b = 42.7e-6
        expected = (parameters['s']+b)/(parameters['s']/30+b/0.3)
        self.assertAlmostEqual(k, expected)

    def test_yang_surface_branch_and_threshold(self):
        boundary = {'mode': 'yang', 'q_actual': 10.0, 'n_coils': 4, 'coil_index': 1,
                    'zeta': 1/0.909, 'psi': 1.0, 'xi': 1.2,
                    'channels': {s: {'area': 0.01, 'diameter': 0.1, 'length': 1.0}
                                 for s in ('inner', 'outer', 'bottom', 'top')}}
        gas = lambda T: (0.1, 1e-5, 0.3, 14000)
        low = surface_coefficients(799, boundary, gas)
        high = surface_coefficients(800, boundary, gas)
        for surface in low:
            self.assertAlmostEqual(high[surface]/low[surface], 1.2)
        boundary['q_actual'] = 0.0001
        with self.assertRaises(ValueError):
            surface_coefficients(800, boundary, gas)

    def test_snapshot_times_and_gas_transition(self):
        settings, grid, material, gas = small_case()
        settings.update(gas_temperatures=[700, 1000], snapshot_interval=137.0)
        result = solve(grid, material, gas, settings)
        self.assertTrue(np.any(np.isclose(result['time'], 200)))  # Tg = 800 K.
        np.testing.assert_allclose(result['field_times'], [0, 137, 274, 411, 548, 600])
        self.assertEqual(result['time'][-1], 600)

    def test_same_material_same_result(self):
        settings, grid, material, gas = small_case()
        a = solve(grid, material, gas, settings)
        b = solve(grid, material, gas, copy.deepcopy(settings))
        np.testing.assert_array_equal(a['fields'], b['fields'])

    def test_target_interpolation_first_crossing(self):
        self.assertEqual(target_time([0, 10, 20], [300, 500, 700], 600), 15)
        self.assertEqual(target_time([0, 10], [700, 800], 600), 0)
        self.assertIsNone(target_time([0, 10], [300, 500], 600))
        self.assertEqual(target_time([0, 10, 20, 30], [300, 700, 400, 800], 600), 7.5)

    def test_unreached_target_has_no_fabricated_difference(self):
        rows = compare_times([{'time': [0, 10], 'cold': [300, 400]},
                              {'time': [0, 10], 'cold': [300, 1000]}], [580])
        self.assertIsNone(rows[0]['time_A_s'])
        self.assertIsNone(rows[0]['difference_A_minus_B_s'])

    def test_bad_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            target_time([0, 0], [300, 400], 350)
        with self.assertRaises(ValueError):
            material_from_table([[250, 7800, 500, 25], [1500, 7700, 500, 25]], 'test', 'test')
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'bad.csv'
            path.write_text('300,7800,500,25\n290,7800,500,25\n')
            with self.assertRaises(ValueError):
                read_table(path, 4)
        settings, grid, material, gas = small_case()
        settings['gas_times'] = [0, 1]
        with self.assertRaises(ValueError):
            solve(grid, material, gas, settings)

    def test_table_mode_with_original_properties_module(self):
        # Isolated synthetic fixtures: never create industrial-looking repo data.
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for path in Path(__file__).resolve().parent.glob('*.py'):
                shutil.copy2(path, root/path.name)
            data = root/'data'
            data.mkdir()
            hydrogen = [[250, 0.10, 1.5e-5, 0.30, 14000], [1500, 0.02, 2e-5, 0.40, 14500]]
            np.savetxt(data/'hydrogen.csv', hydrogen, delimiter=',')
            for name, k in [('steel.csv', 25), ('steel_low.csv', 25), ('steel_high.csv', 45)]:
                np.savetxt(data/name, [[250, 7850, 500, k], [1500, 7850, 800, k]], delimiter=',')
            script = '''
from config import get_config
from main import run_comparison
from material_data import load_models
import numpy as np
s = get_config('tables')
s.update(duration=100, snapshot_interval=100, max_dt=10)
s['geometry'].update(n_r=3, n_z=3)
s['materials'][0]['table'] = 'data/steel.csv'
m, gas = load_models(s)
import properties
np.testing.assert_allclose(m[0]['properties'](500), properties.steel_properties(500))
r, rows = run_comparison(s)
assert len(r) == 2 and r[0]['time'][-1] == 100
assert all(row['time_A_s'] is None for row in rows)
'''
            completed = subprocess.run([sys.executable, '-c', script], cwd=root,
                                       capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == '__main__':
    unittest.main()
