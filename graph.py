"""Matplotlib figures with common axes, scales and explicit data status."""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_temperature_curves(results, targets, output):
    fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
    for result, color in zip(results, ('tab:blue', 'tab:orange')):
        hours = result['time']/3600
        ax.plot(hours, result['cold']-273.15, color=color, label=result['name']+' / cold')
        ax.plot(hours, result['hot']-273.15, color=color, linestyle='--', alpha=0.6,
                label=result['name']+' / hot')
    for target in targets:
        ax.axhline(target, color='gray', linewidth=0.7)
        ax.text(0.99, target, f'{target:g} C', transform=ax.get_yaxis_transform(),
                ha='right', va='bottom', fontsize=8)
    ax.set(xlabel='Time [h]', ylabel='Steel temperature [C]', title=results[0]['status'])
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc='lower right')
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_thermal_maps(results, output):
    reference_times = results[0]['field_times']
    indices = np.unique([min(8, len(reference_times)-1), len(reference_times)//2,
                         len(reference_times)-1])
    fig, axes = plt.subplots(2, len(indices), figsize=(4*len(indices)+1, 7),
                              squeeze=False, layout='constrained', sharex=True, sharey=True)
    lower = min(float(r['fields'].min()) for r in results)-273.15
    upper = max(float(r['fields'].max()) for r in results)-273.15
    if upper <= lower:
        upper = lower+1
    for row, result in enumerate(results):
        for col, index in enumerate(indices):
            time = reference_times[index]
            nearest = int(np.argmin(np.abs(result['field_times']-time)))
            ax = axes[row, col]
            field = result['fields'][nearest]-273.15
            mesh = result['mesh']
            colors = ax.pcolormesh(mesh['r'], mesh['z'], field, shading='gouraud',
                                   cmap='inferno', vmin=lower, vmax=upper)
            cold = np.unravel_index(np.argmin(field), field.shape)
            ax.plot(mesh['r'][cold[1]], mesh['z'][cold[0]], 'cx', markersize=7)
            ax.set(title=f"{result['name']}\nt = {time/3600:.1f} h", xlabel='Radius r [m]')
            if col == 0:
                ax.set_ylabel('Axial position z [m]')
            ax.set_aspect('equal')
    fig.colorbar(colors, ax=axes.ravel().tolist(), label='Temperature [C]', shrink=0.8)
    fig.suptitle(results[0]['status']+'\nSame times and color scale; x marks a coldest node', fontsize=11)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_target_times(results, rows, output):
    fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
    x = np.arange(len(rows))
    for offset, key, result in zip((-0.18, 0.18), ('time_A_s', 'time_B_s'), results):
        values = [np.nan if row[key] is None else row[key]/3600 for row in rows]
        ax.bar(x+offset, values, width=0.36, label=result['name'])
        for i, value in enumerate(values):
            if np.isnan(value):
                ax.text(i+offset, 0, 'Not reached', rotation=90, ha='center', va='bottom', fontsize=8)
    ax.set_xticks(x, [f"{row['target_C']:g}" for row in rows])
    ax.set(xlabel='Target temperature [C]', ylabel='Cold-point first crossing [h]',
           title=results[0]['status'])
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.2)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def save_figures(results, rows, targets, output_directory):
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    plot_temperature_curves(results, targets, output/'temperature_curves.png')
    plot_thermal_maps(results, output/'thermal_maps.png')
    plot_target_times(results, rows, output/'target_times.png')
