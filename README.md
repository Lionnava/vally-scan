# VALLY-Scan v1.4

**[📄 Lea nuestro White Paper para una visión completa del proyecto.](WhitePaper_VALLY.pdf)**

**[🔬 Lea nuestro Manuscrito Científico (Pre-print).](Manuscrito.pdf)**

---

## Manual de Procedimientos: VALLY-Scan v1.4

**Autor:** Lionell E. Nava Ramos  
**Versión:** 1.4 (Prototipo para el Premio Nacional de Inventiva Tecnológica 2025)  
**Fecha:** Octubre 2024

### 1. Introducción y Propósito

El software **VALLY-Scan v1.4** es la primera implementación del marco computacional desarrollado bajo el **Proyecto V.A.L.L.Y. (Vibrational Analysis for Ligand Likelihood Yielding)**.

**Propósito:** Proporcionar a los investigadores una herramienta de línea de comandos para realizar un análisis rápido y preliminar de la dinámica de complejos proteína-ligando. El sistema calcula la flexibilidad intrínseca de la proteína mediante un Modelo de Red Elástica (ENM) y utiliza esta información para predecir una afinidad de enlace plausible y los residuos de mayor impacto en la interacción, todo en cuestión de segundos y en hardware de escritorio estándar.

### 2. Requisitos del Sistema

*   **Sistema Operativo:** Windows, macOS o Linux.
*   **Software:** Distribución de Python a través de Anaconda o Miniconda.
*   **Hardware Mínimo:** Procesador Core i3 o equivalente, 8 GB de memoria RAM.
*   **Conexión a Internet:** Requerida únicamente para la instalación inicial de las dependencias.

### 3. Procedimiento de Implementación (Instalación desde Cero)

#### Paso 3.1: Clonar o Descargar el Repositorio
La forma más fácil es descargar el proyecto como un archivo ZIP. En la página principal de este repositorio, haz clic en el botón verde `<> Code` y luego en `Download ZIP`. Descomprime el archivo en tu computadora.

#### Paso 3.2: Creación del Entorno Virtual Conda
Abra la terminal de comandos (Anaconda Prompt en Windows).
Ejecute el siguiente comando para crear un entorno aislado llamado `marco-bio-ia` con Python 3.9:
```bash
conda create --name marco-bio-ia python=3.9

Active el entorno recién creado:
code
Bash
conda activate marco-bio-ia
Paso 3.3: Instalación de Dependencias
Con el entorno activado, proceda a instalar las librerías de Python requeridas.
Nota para Windows: Si se presentan errores de certificados, ejecute set CURL_CA_BUNDLE= antes de los siguientes comandos.
code
Bash
pip install numpy prody reportlab matplotlib
4. Puesta en Marcha (Ejecución del Análisis)
Una vez implementado, el sistema se ejecuta a través de la terminal de comandos.
Paso 4.1: Iniciar Sesión de Trabajo
Abra la terminal de comandos (Anaconda Prompt).
Active el entorno de trabajo:
code
Bash
conda activate marco-bio-ia
Nota para Windows: Si aplica, ejecute set CURL_CA_BUNDLE=.
Navegue al directorio donde descomprimió el proyecto (la carpeta vally-scan):
code
Bash
cd ruta/a/la/carpeta/vally-scan
Paso 4.2: Ejecución del Análisis
Para ejecutar un análisis, utilice el siguiente comando, reemplazando [pdb_id] por el código de 4 letras del archivo PDB a analizar (ej. 2fom o 6lu7).
code
Bash
python main.py --pdb_id [pdb_id]
Para generar adicionalmente un reporte en PDF, añada el argumento --reporte.
code
Bash
python main.py --pdb_id [pdb_id] --reporte
El reporte será guardado en una nueva carpeta llamada reportes.
Ejemplos de Uso:
code
Bash
python main.py --pdb_id 2fom
code
Bash
python main.py --pdb_id 6lu7 --reporte
Paso 4.3: Visualización de Ayuda
Para ver todas las opciones y la descripción del programa, ejecute:
code
Bash
python main.py --help
5. Licencia
Este proyecto está distribuido bajo la licencia MIT. Consulte el archivo LICENSE para más detalles.
