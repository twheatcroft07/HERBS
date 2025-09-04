# HERBS
A Python-based GUI for Histological E-data Registration in Brain Space

HERBS is an open source, extensible, intuitive and interactive software platform for image visualisation and image registration. Image registration is the process of identifying a spatial transformation that maps images to a template such that corresponding anatomical structures are optimally aligned — in other words, a voxel-wise ‘correspondence’ is established between the images and template.

HERBS has been tested on Windows 10, macOS (Big Sur – Monterey), and Linux (Kubuntu 18.04, Ubuntu 22.04 LTS). It runs reliably in Python 3.9.0–3.9.2 with PyQt5 ≥ 5.14.2 as a GUI framework.

---

## 🔧 TW Recommended Changes

<span style="color:darkred; font-weight:bold;">Based on known issues with Python ≥3.10 (see [#71](https://github.com/JingyiGF/HERBS/issues/71), [#68](https://github.com/JingyiGF/HERBS/issues/68), [#39](https://github.com/JingyiGF/HERBS/issues/39), [#25](https://github.com/JingyiGF/HERBS/issues/25)), we recommend using **Python 3.9.2**. This avoids crashes during atlas/mesh download and processing.</span>

---

## Install (TW Recommended)

We strongly recommend using the provided `environment.yml` file.  
This pins Python 3.9.2 and dependency versions that are known to work with atlas/mesh downloads.

```bash
# 1. Clone the TW fork
git clone https://github.com/twheatcroft07/HERBS.git
cd HERBS

# 2. Create the environment
conda env create -f environment.yml

# 3. Activate the environment
conda activate herbs39

# 4. Install HERBS in editable mode
pip install -e .
```

That’s it — you now have HERBS installed in a reproducible, working environment.

---

## Usage

Activate your environment first:
```bash
conda activate herbs39
```

Then run HERBS:
```python
import herbs
herbs.run_herbs()
```

A GUI window will open. Users can download or manually load atlases, then upload images for further processing.

<img src="./herbs/herbs.png" width="800px"></img>

For more information, please read the HERBS CookBook (ongoing) or check the Tutorial folder.

---

## Atlas Storage and Download

- **Do not** store atlases inside the HERBS folder itself.
- Store each atlas version in its own separate folder.

⚠️ Known issue: the built-in atlas downloader may hang or fail on Python ≥3.10 (see also [#71](https://github.com/JingyiGF/HERBS/issues/71), [#68](https://github.com/JingyiGF/HERBS/issues/68)).

✅ With **Python 3.9.2**, we do *not* expect this issue — atlas downloads and processing work reliably in our testing.

- **Manual atlas setup (alternative or backup method):**
  1. Create a folder outside the HERBS repo (e.g. `~/HERBS_ATLASES/allen_ccf_10um`).
  2. Place the Allen atlas files there:
     - `annotation_10.nrrd`
     - `average_template_10.nrrd`
     - (plus optional but recommended: `atlas_labels.pkl`, `atlas_axis_info.pkl`, `atlas_meshdata.pkl`, etc.)
  3. In HERBS, select this folder when prompted, then click *Process*.

---

## Requirements & Dependencies

- HERBS requires a 64-bit operating system and 64-bit Python.
- 3D visualisation in HERBS depends on OpenGL. If you encounter OpenGL errors:
  - macOS: edit `OpenGL/platform/ctypesloader.py` to set the path to `/System/Library/Frameworks/OpenGL.framework/OpenGL`.
  - Linux/Windows: install OpenGL from [https://www.opengl.org](https://www.opengl.org).

## Dependencies (TW Recommended / Tested)

HERBS runs best in a conda environment built from the provided [`environment.yml`](./environment.yml).  
This YAML pins the tricky packages (NumPy, pandas, Numba, PyQt) to compatible versions and ensures all required plotting/GUI libraries are included.

- **Python**: 3.9.2  
- **NumPy**: 1.23.5  _(avoid 2.x; pyqtgraph paths use `np.product`)_  
- **pandas**: 1.5.3  _(built against NumPy 1.23 ABI)_  
- **Numba**: 0.56.4  _(compatible with NumPy 1.23)_  
- **PyQt / PyQt5**: Qt 5.x / PyQt5 ≥ 5.14.2  
- **pyqtgraph**: 0.12.3  
- **matplotlib**  
- **scikit-image**  
- **scipy**  
- **plotly**  
- **dash**

---

**Notes**  
- **NumPy 2.x** is not supported (removes `np.product`, breaks rotated slice views).  
- If you generated atlas pickle files under a different NumPy major version, re-run **Process** in HERBS to rebuild them.  
- Installing via the `environment.yml` is strongly recommended to ensure consistent versions across machines.


If dependency conflicts arise, the easiest fix is to create a fresh environment and reinstall HERBS.

---

## Issues & Discussions

- Report issues: [https://github.com/JingyiGF/HERBS/issues](https://github.com/JingyiGF/HERBS/issues)  
- Start discussions: [https://github.com/JingyiGF/HERBS/discussions](https://github.com/JingyiGF/HERBS/discussions)

When reporting, please include a clear description, error messages, or screenshots.

---

## Development Notes

HERBS is under active development. Check for updates frequently, especially before starting a new project. Fork users can push changes to their repo while still being able to pull upstream updates if needed.

---

Hope this tool makes your research life more tasty :-)

