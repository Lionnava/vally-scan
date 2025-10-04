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
