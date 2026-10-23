# Recocido: comparación axisimétrica de dos materiales

**El modo predeterminado es ilustrativo. No contiene propiedades ni dimensiones
verificadas de Ternium y sus resultados no son tiempos industriales validados.**

Inspección del repositorio: commit `ecc05837972e36c8e094a35826ee20528c07979c`.
Se conservaron sin cambios `properties.py`, `flow.py`, `convection.py`, `mesh.py`
y `radial_conductivity.py`. En ese commit **no existía la carpeta `data/`**.

## Ejecutar en VS Code (Windows)

Instalar Python 3.10 o posterior y la extensión Python. Abrir esta carpeta en
VS Code y abrir una terminal PowerShell en la carpeta que contiene `main.py`:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m unittest -v test_model
.\.venv\Scripts\python.exe main.py --mode demo --output results_demo
```

No hace falta activar el entorno. Para usar el botón de ejecución de VS Code,
seleccionar `.venv\Scripts\python.exe` en **Python: Select Interpreter**.
En macOS/Linux, usar `python3 -m venv .venv` y `.venv/bin/python` en lugar de la
ruta de Windows. El programa guarda las gráficas como PNG; no abre ventanas.

## Ejecutar en Google Colab

Usar el ZIP completo entregado (`Recocido_pruebas_completo.zip`), porque los
archivos nuevos todavía no están publicados en GitHub. En una primera celda:

```python
from google.colab import files
import zipfile
files.upload()  # Seleccionar Recocido_pruebas_completo.zip.
with zipfile.ZipFile('Recocido_pruebas_completo.zip') as bundle:
    bundle.extractall('/content')
```

En otra celda:

```python
%cd /content/Recocido_pruebas
%pip install -r requirements.txt
!python -m unittest -v test_model
!python main.py --mode demo --output results_demo
```

Para visualizar y descargar resultados:

```python
from IPython.display import Image, display
display(Image('results_demo/temperature_curves.png'))
display(Image('results_demo/thermal_maps.png'))
import shutil
shutil.make_archive('/content/resultados_recocido', 'zip', 'results_demo')
files.download('/content/resultados_recocido.zip')
```

## Archivos y conexiones

| Archivo nuevo | Responsabilidad y dependencias |
|---|---|
| `config.py` | Geometría, ciclo del gas, materiales, contacto, fronteras y ajustes. |
| `material_data.py` | Propiedades demo o CSV; adapta `properties.py` y proporciona entalpía e inversa. |
| `surface_transfer.py` | Coeficientes prescritos o Yang (5)-(8); usa `flow.surface_flows` y `convection.convection_heat_flux`. |
| `solver.py` | Usa `mesh.stable_time_step`, `radial_conductivity.effective_radial_conductivity` y los adaptadores. |
| `analysis.py` | Primer cruce del punto frío, diferencia A-B y balance energético. |
| `graph.py` | Curvas, mapas con escala común y tiempos objetivo mediante Matplotlib. |
| `main.py` | Construye la malla con `mesh.build_mesh`, ejecuta ambos materiales y exporta. |
| `test_model.py` | Pruebas independientes, físicas y de integración con unittest. |
| `verify_convergence.py` | Repite la comparación con menor paso y mallas más finas. |

Para incorporar el desarrollo al repositorio del equipo, copiar los archivos
nuevos del ZIP `Recocido_pruebas_solo_nuevos.zip`; no reemplaza los módulos
originales. En ese ZIP la documentación se llama `README_RUN.md` para preservar
también el README del equipo. El ZIP completo puede ejecutarse en otra carpeta.

## Salidas y criterio de tiempo

`results_demo/` contiene tres PNG, `target_times.csv`, dos historias CSV, dos
archivos NPZ con campos completos y `run_metadata.json` con parámetros y diagnósticos.
Las historias usan segundos, kelvin, metros y joules para el rollo completo.
`target_times.csv` usa grados Celsius para los objetivos y segundos para tiempos.
Una celda vacía indica que no se alcanzó el objetivo; no significa tiempo cero.

Se calcula el primer instante en que **todos los nodos** alcanzan el objetivo:
`min(T) >= T_obj`. El cruce se interpola entre pasos consecutivos, no entre las
instantáneas gráficas. Las posiciones del punto frío y caliente se guardan en
cada paso. Con varios nodos empatados se registra uno de ellos.

## Modelo y unidades

- Sección `(n_z, n_r)` con nodos sobre las superficies: fila 0 inferior y columna
  0 interior. Áreas y volúmenes vienen directamente de `mesh.py`, por radián.
- Balance conservativo axisimétrico: conducción radial y axial, sin generación.
  Las caras internas intercambian potencias iguales y opuestas. Los nodos de
  esquina reciben intercambio por ambas superficies.
- Único modelo radial: función existente de Yang (9)-(17), que usa `T_m**3`
  coherentemente con la conductancia radiativa de (11). No se sustituyó.
- `k_z` es la conductividad sólida; `k_r` es la efectiva calculada. Se evalúa la
  conductividad del gas interlaminar a la temperatura local del acero. `T_g`
  selecciona la rama de 800 K del módulo existente.
- Euler explícito en entalpía volumétrica, con inversión exacta para `cp(T)`
  lineal por tramos. Densidad constante, sin expansión. Conductancias armónicas
  en las caras. Se recalculan coeficientes y límite de paso en cada iteración.
- Se usa el mínimo tabulado de `rho*cp` en el límite de estabilidad para cubrir
  también segmentos donde `cp` decrece. Esto puede ser conservador.
- Fronteras Robin `q''=h_eff*(T_g-T_s)`, positivas hacia el rollo. `h=0` representa
  una superficie adiabática. El gas es uniforme y prescrito, sin balance de gas.
- Por defecto, `h_eff=80 W/(m2 K)` en las cuatro superficies es **provisional**.
  No se añade otra radiación exterior a ese coeficiente. La radiación interlaminar
  permanece dentro del modelo radial y es un intercambio distinto.
- El modo opcional `boundary['mode']='yang'` requiere caudal circulante real,
  sección de paso, diámetro hidráulico, longitud, `psi` y `xi` documentados.
  No basta con los diámetros del rollo. Se rechaza Re < 10000; superar ese valor
  no demuestra por sí solo que la correlación esté validada para ese canal.
- Presión del gas en Pa; presión de contacto y dureza en MPa; rugosidad y espesores
  en m. Temperaturas internas en K. La transición de Yang es 800 K, no 800 C.

## Pasar a datos reales y verificar

Consultar `data/README.md`. Completar las tablas con fuentes, nombres de grados,
geometría real, ciclo, parámetros de contacto y fronteras en `config.py`.
`--mode tables` nunca genera tablas ni sustituye datos ausentes silenciosamente:

```powershell
.\.venv\Scripts\python.exe main.py --mode tables --output results_tables
.\.venv\Scripts\python.exe verify_convergence.py
```

El segundo comando estudia **el caso demo** y guarda `results_demo/convergence.json`.
Para estudiar otros parámetros, cambiar explícitamente su llamada a `get_config`.
También se puede ejecutar `main.py --nr 25 --nz 33 --max-dt 10 --output results_fine`.
La conservación y la convergencia verifican el programa; **no validan el horno**.

La demostración mantiene iguales geometría, `rho`, `cp`, fronteras y contacto,
y cambia únicamente `k_s`: 25 frente a 45 W/(m K). El gas demo tiene propiedades
simplificadas explícitas, no una tabla atribuida a CoolProp. Ambos casos pasan por
la misma función radial. Todos estos parámetros requieren revisión industrial.
Comparar siempre en la misma posición de la carga. No asumir que el orden de las
difusividades persiste con tablas reales ni que el tiempo depende solo de alfa.

Referencia principal: Yang et al. (2025),
https://doi.org/10.1038/s41598-025-97422-4.
Fechas del proyecto: exposición 23/10/2026; reporte 06/11/2026.
