# Verificación del desarrollo - 22 de octubre de 2026

## Alcance inspeccionado

Repositorio: https://github.com/Acsl15/Recocido_pruebas.git

Commit: `ecc05837972e36c8e094a35826ee20528c07979c`.

| Archivo original | Estado observado y tratamiento |
|---|---|
| properties.py | Interpoladores definidos, pero faltaban ambos CSV. Se conserva; se usa al cargar tablas verificables. |
| flow.py | Reproduce (5)-(6) para cuatro rollos con zeta=1/0.909. Se conserva. |
| convection.py | Calcula flujo/potencia a partir de h; no calcula h. Se conserva y se llama desde surface_transfer.py. |
| mesh.py | Malla por radián, arreglos (n_z,n_r), áreas y límite explícito utilizables. Se conserva. |
| radial_conductivity.py | Modelo Yang con mecanismos en paralelo y conductancia radiativa proporcional a T^3. Se conserva íntegro. |
| data/ | No figuraba en el árbol del commit. Se añade solamente documentación; no se fabrican tablas industriales. |
| solver, analysis, graph, config, main | No existían; se añadieron. |

No se publicaron commits ni cambios en GitHub. Los cinco archivos Python originales
permanecen sin diferencias respecto al repositorio inspeccionado.

## Pruebas ejecutadas

19 pruebas independientes y de integración aprobadas, incluyendo:

- Volumen y cuatro áreas del cilindro frente a sus fórmulas geométricas exactas.
- Distribución de flujo frente a la tabla de Yang y conservación del caudal.
- Cancelación de los intercambios internos y potencia superficial del cilindro.
- Entalpía contra la integral analítica de cp lineal; inversión con cp decreciente.
- Equilibrio uniforme; límites de temperatura durante calentamiento y enfriamiento.
- Energía calculada independientemente con cp constante y con cp variable.
- Decaimiento de un modo coseno axial aislado contra la solución analítica, con
  reducción de error al refinar la malla.
- Límite de resistencias en serie del modelo radial original.
- Rama convectiva Yang, cambio a 800 K y rechazo de un régimen no turbulento.
- Instantes de salida y transición térmica; repetibilidad para un mismo material.
- Primer cruce interpolado, objetivo inicialmente alcanzado y objetivo no alcanzado.
- Rechazo de entradas inválidas; integración con properties.py original usando
  CSV sintéticos en una carpeta temporal independiente, que no se distribuyen
  como datos industriales.

Las comprobaciones reproducen también una ejecución completa con gráficas.
Las tres imágenes se revisaron visualmente: ejes legibles, etiquetas y escala
común de temperatura en los mapas.

## Comparación exclusivamente ilustrativa

Diámetros 0.60/1.60 m; ancho 1.00 m; temperatura inicial 25 C; gas de 25 a 800 C
en cuatro horas y después constante; 24 horas simuladas; h_eff=80 W/(m2 K) en
las cuatro superficies. Rho=7850 kg/m3 y cp=600 J/(kg K) iguales en ambos casos.
Solo cambia la conductividad sólida: A=25 y B=45 W/(m K). La conductividad radial
se calcula en ambos mediante el módulo Yang existente. Los demás parámetros
provisionales están en config.py y run_metadata.json.

| Objetivo del punto frío | Caso A [h] | Caso B [h] | A-B [h] |
|---|---:|---:|---:|
| 580 C | 9.300 | 7.975 | 1.324 |
| 680 C | 12.052 | 10.312 | 1.739 |
| 780 C | 20.167 | 17.210 | 2.957 |

Tiempo por interpolación del primer cruce del mínimo nodal. No incluye permanencia
metalúrgica ni enfriamiento. Estas cifras no representan ningún grado de Ternium.

## Estabilidad, energía y sensibilidad numérica

- Malla base: 13 nodos radiales x 17 axiales. Paso máximo 20 s, limitado además
  por el balance de cada volumen con factor de seguridad 0.85.
- 4321 pasos por material. Ningún paso excede su límite calculado.
- Residuo energético relativo máximo de la corrida base: aproximadamente
  2.04e-15, normalizado por la máxima energía neta acumulada desde las fronteras.
  Este valor verifica el balance discreto, no la precisión industrial.
- Reducir el máximo de 20 a 10 s cambia los tiempos objetivo como máximo 0.02310 %.
- A igual máximo de 10 s, pasar de 13x17 a 25x33 cambia los tiempos como máximo
  0.01892 %; pasar de 25x33 a 37x49 cambia como máximo 0.003503 %.
- Los detalles de los ocho cálculos están en results_demo/convergence.json.
  Son diferencias observadas en este caso, no una cota universal del error.

## Pendientes para una comparación industrial válida

1. Identificar dos grados reales y justificar su orden de difusividad en el rango
   completo. Agregar rho, cp(T), k(T) y procedencia; verificar posibles cruces.
2. Confirmar dimensiones, espesor de lámina, presión de contacto, rugosidad,
   parámetros de óxido/emulsión y aplicabilidad de Yang a los rollos reales.
3. Incorporar ciclo de gas y coeficientes medidos o, para (7)-(8), geometría de
   canales, caudal circulante real y factores psi/xi documentados. No duplicar
   radiación exterior si ya está en h_eff.
4. Mantener la misma geometría, posición de carga y condiciones externas para
   ambos materiales; distinguir punto frío de temperatura de gas y objetivo.
5. Repetir convergencia con los datos reales y contrastar con mediciones. El
   estudio actual verifica funcionamiento numérico, no valida un horno real.

## Límites deliberados

Funciones sencillas y NumPy, sin clases de aplicación ni dependencias de solver.
Un único rollo por simulación, medio continuo axisimétrico anisótropo, densidad
constante, gas uniforme prescrito, sin cambios geométricos, cinética metalúrgica
ni balance acoplado de gas. Las propiedades fuera del rango de tabla se rechazan.
El código no corrige las fórmulas del equipo ni introduce un modelo radial alterno.
