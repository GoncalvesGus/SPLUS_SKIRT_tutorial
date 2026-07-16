# SPLUS + SKIRT Tutorial

Hands-on tutorial covering the full pipeline from raw IllustrisTNG (TNG50) particle
data through SKIRT radiative transfer to synthetic S-PLUS mock images with realistic
instrumental noise and PSF.

**Author:** Gustavo Fernandes Gonçalves
ORCID: [0009-0006-5887-6621](https://orcid.org/0009-0006-5887-6621)
Lattes: http://lattes.cnpq.br/7416443230735445

<img width="500" height="259" alt="output" src="https://github.com/user-attachments/assets/a9bd17f2-8345-4f02-b195-74ab7030c1a1" />
<img width="1433" height="863" alt="Data_analysis" src="https://github.com/user-attachments/assets/470efecf-e6c8-40f3-80f1-6a929a1ee96a" />



## Prerequisites

Python 3 installed. Before the tutorial, create a virtual environment and install
JupyterLab in it:

```bash
python3 -m venv splus_skirt_env
source splus_skirt_env/bin/activate
pip install jupyterlab
```

You'll need to run `source splus_skirt_env/bin/activate` again each time you open a
new terminal to work on this tutorial.

## Getting started

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/SPLUS_SKIRT_tutorial.git
   cd SPLUS_SKIRT_tutorial
   ```
2. Download the large data files (SKIRT datacubes and TNG50 cutouts) from the
   GitHub Release:
   ```bash
   bash download_data.sh
   ```
3. Start JupyterLab and open the notebooks in order:
   ```bash
   jupyter lab
   ```
   - `Notebook_0__setup_install.ipynb` — checks/installs required Python libraries and DS9
   - `Notebook_1__raw_TNG50_data.ipynb` — inspecting raw TNG50 particle data
   - `Notebook_2__SKIRT_data_SPLUS_filters.ipynb` — SKIRT datacube inspection and S-PLUS filter convolution
   - `Notebook_3__mock_noises_PSF.ipynb` — adding instrumental realism (resampling, noise, PSF)

## Repository structure

```
Notebook_0__setup_install.ipynb
Notebook_1__raw_TNG50_data.ipynb
Notebook_2__SKIRT_data_SPLUS_filters.ipynb
Notebook_3__mock_noises_PSF.ipynb
tutorial_tools.py             # shared helper functions
download_data.sh              # downloads large data from the GitHub Release
test_42.fits                  # small test FITS used in Notebook 0
filter_curves-master/         # S-PLUS filter transmission curves
SKIRT_datacubes/               # SKIRT output datacubes (downloaded, not committed)
TNG50_hdf5_cutouts/            # TNG50 raw + Firefly-ready cutouts (downloaded, not committed)
```

`fits_noisefree/` and `fits_mock_images/` are created automatically as you run the
notebooks.

## License

MIT
