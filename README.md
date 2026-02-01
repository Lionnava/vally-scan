V.A.L.L.Y. Project: Universal Computational Framework (v1.7)
Redefining Agility in Drug Discovery and Materials Science

🔬 Descripción General
V.A.L.L.Y. (Vibrational Analysis for Ligand Likelihood Yielding) es un marco computacional de ultra-alta velocidad diseñado para identificar sitios dinámicos críticos en estructuras complejas. Mediante la integración de Modelos de Red Anisotrópica (ANM) y una lógica de Validación de Triple Factor, el framework conecta la física teórica con la cristalografía experimental.

🚀 Innovaciones de la v1.7 (Universal Edition)
Esta versión marca el salto de una herramienta bioinformática a un motor de física universal:Validación de Triple Factor (Sistema R2):Factor 1 (Física): Cálculo de la dinámica intrínseca mediante la Matriz Hessiana.Factor 2 (Heurística Informada): Discriminación automática entre sitios catalíticos y alostéricos por exclusión geométrica.Factor 3 (Realidad Experimental): Correlación en tiempo real ($r$) con B-factors de cristalografía.Núcleo Universal: Soporte multimodal para Proteasas Virales (SARS-CoV-2, Dengue) y Estructuras Minerales/Cristalinas.Rendimiento Extremo: Análisis completados en segundos en hardware estándar (Core i3, 8GB RAM), democratizando el acceso a la biofísica avanzada.

🛠️ Estructura del Repositorio
vally_scan_v1_7_universal.py: El motor principal con el sistema de triple validación.

main.py: Versión estable de producción.

data/: Archivos PDB y CIF validados (6LU7, 2FOM).

plots/: Gráficas de validación que demuestran la correlación con cristalografía.

📊 Resultados ValidadosEl framework ha sido validado contra objetivos virales de alta prioridad:| Objetivo | PDB ID | Correlación Pearson ($r$) | Descubrimiento Clave || :--- | :--- | :--- | :--- || SARS-CoV-2 Mpro | 6LU7 | 0.67 | Cluster Alostérico detectado (277-279) || Dengue Protease | 2FOM | 0.43 | Región de alta flexibilidad (43-45) |

🎓 Autor e Institución
Ing. Lionell E. Nava Ramos
Investigador en el Centro de Investigación en Informática (CII)
Universidad Politécnica Territorial de Maracaibo (UPTMA)
Ponente aceptado en la APS Global Physics Summit 2026, Denver, CO.
