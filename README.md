# Morelia

Morelia is a wxPython desktop application that collects molecular-modelling,
trajectory-analysis, conversion, visualisation, docking, and plotting utilities
behind one menu-driven interface. It is primarily designed for Linux and
integrates with several external scientific applications.

The main program is `bin/morelia.py`.

## Functionality

Morelia provides the following groups of tools:

- **Desmond**
  - trajectory alignment;
  - ligand and protein RMSD analysis;
  - protein RMSF and protein-ligand contact analysis;
  - Desmond minimisation and molecular-dynamics launchers;
  - Simulation Interaction Diagram, Simulation Quality Analysis, and trajectory
    clustering launchers.
- **GROMACS**
  - `.trr` to `.xtc` trajectory compression with `gmx trjconv`;
  - a Desmond-to-gmx_MMPBSA preparation and calculation workflow.
- **Python**
  - measurement of ligand-atom distances from receptor residue centroids.
- **VMD**
  - Desmond `.cms` to GROMACS `.gro` conversion;
  - Desmond `.dtr` trajectory to GROMACS `.trr` conversion.
- **YASARA**
  - docking, rescoring, ligand preparation, virtual-screening helpers;
  - hydrogen-bond, hydrophobic, ionic, pi-cation, and pi-pi interaction
    analysis;
  - measurement of distances between a receptor residue atom and a ligand
    atom across docking poses (`yasara_distance_app.py`); pose sorting
    optionally appends these distances as the last column of the sorted
    results table;
  - PDB/membrane alignment and YOb conversion/inspection tools.
- **Results and visualisation**
  - Grace plotting and graph viewing;
  - table viewing and heat plots;
  - PDB, MOL2, and CMS viewing through PyMOL;
  - PDB sequence and YOb result inspection.

Many tools operate on files in the current working directory and create result
files there. Select the intended directory with **File → Working Directory**
before starting a calculation.

## Requirements

### Operating system and Python

- Linux is recommended. A number of legacy utilities invoke Unix shell commands
  and Linux desktop applications directly.
- Python 3.10 or newer is recommended.
- wxPython is required for the graphical interface.

The directly imported third-party Python packages are listed in
`requirements.txt`:

- Biopython
- Matplotlib
- MDTraj
- NetworkX
- NumPy
- pandas
- ParmEd
- wxPython

The YASARA example also imports `yasaramodule`. That module is supplied by
YASARA and is not installed from PyPI.

### External applications

Only install the applications required for the Morelia features you intend to
use.

| Application | Used for | Expected command or configuration |
| --- | --- | --- |
| Maestro and Desmond | Desmond setup, simulation, analysis, and clustering tools | `maestro` must be available; set the `SCHRODINGER` environment variable so `$SCHRODINGER/run` works |
| GROMACS | `.trr` to `.xtc` compression and MM-PBSA preparation/calculation | `gmx` in `PATH` |
| VMD | `.cms` to `.gro` and `.dtr` to `.trr` conversion | `vmd` in `PATH`, including plugins capable of reading Maestro/Desmond and DTR data |
| YASARA | Docking, interaction analysis, YOb tools, and macros | `yasara` in `PATH`; `yasaramodule` is needed for the Python example |
| Grace-compatible plotting program | Interactive and batch plotting | Configure `ggrace`, `xmgrace`, `gracegtk`, `gracebat`, or equivalent commands |
| PyMOL | Optional PDB, MOL2, and CMS visualisation | `pymol` in `PATH` |
| IPython | Optional interactive console launcher | `ipython3` in `PATH` |

Maestro/Desmond and YASARA are separately distributed and licensed products.
Install them according to their vendors' documentation and ensure that your
licence permits the intended use.

#### Additional MM-PBSA requirements

The **GROMACS → Desmond MM-PBSA** mini-application requires a suitable Conda
environment or an equivalent active Python environment containing:

- NumPy, MDTraj, NetworkX, and ParmEd;
- GROMACS (`gmx`);
- `gmx_MMPBSA`;
- AmberTools executables `parmchk2`, `tleap`, and `antechamber`;
- an MPI launcher such as `mpirun`.

The mini-application includes **Check environment** under its Tools menu. Use
that check before starting a preparation or calculation.

## Installation

1. Enter the project directory:

   ```bash
   cd /path/to/morelia
   ```

2. Create and activate a virtual environment:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   ```

3. If wxPython cannot be installed from PyPI, use your operating system's
   package instead. On Debian or Ubuntu:

   ```bash
   sudo apt install python3-wxgtk4.0
   ```

   One practical option is to recreate the virtual environment with access to
   system packages and install the remaining libraries there:

   ```bash
   python3 -m venv --system-site-packages .venv
   source .venv/bin/activate
   python -m pip install biopython matplotlib mdtraj networkx numpy pandas ParmEd
   ```

4. Install the external applications needed for your workflows and make their
   command-line launchers available in `PATH`.

## Configuration

Edit `bin/preferences.py` before first use:

- set `install_dir` to the absolute path of the `bin/` directory;
- set `doc_dir` and `example_dir` to the desired documentation and example
  directories;
- set the commands used to launch Grace, Maestro, VMD, YASARA, IPython, and a
  terminal.

For example:

```python
install_dir = "/path/to/morelia/bin/"
example_dir = "/path/to/morelia/Examples/"

cmdGG = "xmgrace"
cmdGB = "gracebat"
cmdMaestro = "maestro"
cmdVMD = "vmd"
cmdYASARA = "yasara"
```

Grace needs both an interactive command (`cmdGG`) and a batch command (`cmdGB`).
The defaults use `ggrace` and `ggrace -batch`; alternatives shown in
`preferences.py` include `xmgrace`, `gracegtk`, and `gracebat`.

For Maestro/Desmond, initialise the Schrödinger environment before launching
Morelia. A typical shell configuration resembles:

```bash
export SCHRODINGER=/path/to/schrodinger
export PATH="$SCHRODINGER:$PATH"
```

Use the path supplied by your local Schrödinger installation.

## Running Morelia

From the project root:

```bash
source .venv/bin/activate
python bin/morelia.py
```

The **Applications** menu can launch the configured external programs. The
**Preferences** item opens an editor for `bin/preferences.py`.

## Project layout

```text
morelia/
├── bin/                  Main application, mini-applications, and YASARA macros
├── Examples/             Example structures, scripts, and data
├── requirements.txt      Direct Python dependencies
└── README.md
```

## Troubleshooting

- **`ModuleNotFoundError: No module named 'wx'`**  
  Install wxPython for the same Python interpreter used to launch Morelia.
- **An external program does not open**  
  Run its configured command in a terminal and update `bin/preferences.py` if
  necessary.
- **Desmond tools do not start**  
  Check that `SCHRODINGER` is set and that `$SCHRODINGER/run` exists.
- **A GROMACS conversion fails**  
  Check that `gmx` is in `PATH` and that the expected matching structure and
  trajectory files exist.
- **A VMD conversion fails**  
  Check `vmd -dispdev text` and confirm that the installed VMD plugins can read
  the relevant Maestro/Desmond format.
- **A YASARA operation fails**  
  Check the configured YASARA command and confirm that the macros under
  `bin/YASARA/` are accessible.
- **Grace opens but batch plots fail**  
  Verify both `cmdGG` and `cmdGB`; interactive and batch commands may have
  different executable names.

## Licence

The application's About dialog identifies Morelia as non-commercial software
under a Creative Commons BY-NC licence. Add a standalone `LICENSE` file if the
project is to be redistributed.

## Author

Jan Jakubik — `jan.jakubik@fgu.cas.cz`
