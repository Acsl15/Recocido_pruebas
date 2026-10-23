# Datos pendientes

Esta carpeta no contenía tablas en el commit inspeccionado. No se han creado
CSV ficticios con apariencia de datos industriales. `--mode demo` no los necesita.

Para `--mode tables`, agregar CSV numéricos, sin encabezado o con encabezado `#`:

- `hydrogen.csv`: T [K], rho [kg/m3], mu [Pa s], k [W/(m K)], cp [J/(kg K)].
- `steel.csv`: T [K], rho [kg/m3], cp [J/(kg K)], k sólido [W/(m K)]. El módulo
  original lo carga al importar aunque se seleccionen otras tablas.
- `steel_low.csv` y `steel_high.csv`: el mismo formato de cuatro columnas.
  Alternativamente, apuntar una entrada `materials[*]['table']` a `steel.csv`.

Temperaturas estrictamente crecientes, al menos dos filas, propiedades positivas,
sin valores ausentes. Rho debe ser constante: no se modela expansión ni movimiento
de la malla. Las tablas deben cubrir el rango completo de gas y acero, desde la
condición inicial, no solamente 580-780 C. No hay extrapolación.

Documentar grado, composición, fuente, unidades, intervalo y presión del hidrógeno.
El módulo original usa una atmósfera; otra presión requiere actualizar de forma
coherente la tabla y la definición de caudal. No mezclar Pa del gas con MPa de
contacto entre capas. Renombrar los materiales en config.py al cargar datos reales.
