"""
splus_tools.py - auxiliary functions for the SKIRT + SPLUS tutorial
"""

import subprocess
import time

import numpy as np
import astropy.units as u
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.visualization import AsinhStretch


trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz  # numpy >=2.0 removed trapz

# ---------------------------------------------------------------------------
# Notebook 1 — raw TNG50 data
# ---------------------------------------------------------------------------
 
FIELD_DESCRIPTIONS = {
    "Coordinates": "Particle position (ckpc/h, comoving)",
    "Velocities": "Peculiar velocity (km/s * sqrt(a))",
    "Masses": "Particle mass (1e10 Msun/h)",
    "ParticleIDs": "Unique particle identifier",
    "Density": "Local gas density (1e10 Msun/h / (ckpc/h)^3)",
    "InternalEnergy": "Gas specific internal energy (used for temperature)",
    "ElectronAbundance": "Free electron fraction (ne/nH)",
    "StarFormationRate": "Star formation rate of the gas cell (Msun/yr)",
    "GFM_Metallicity": "Mass fraction in metals (Z)",
    "GFM_InitialMass": "Initial stellar particle mass at birth (1e10 Msun/h)",
    "GFM_StellarFormationTime": "Scale factor 'a' at formation (>0: real star; <=0: wind particle)",
}
 
 
def summarize_fields(path, ptype_key):
    """HDF5 particle group -> DataFrame (field, shape, dtype, description)."""
    import h5py
    import pandas as pd
 
    with h5py.File(path, "r") as f:
        group = f[ptype_key]
        rows = [{
            "field": name,
            "shape": dset.shape,
            "dtype": str(dset.dtype),
            "description": FIELD_DESCRIPTIONS.get(name, "-"),
        } for name, dset in group.items()]
    return pd.DataFrame(rows).sort_values("field").reset_index(drop=True)
 
 
def preview_particles(path, ptype_key, n=5):
    """First n particles of each field -> DataFrame (3-vectors split into _x/_y/_z)."""
    import h5py
    import pandas as pd
 
    columns = {}
    with h5py.File(path, "r") as f:
        group = f[ptype_key]
        for name, dset in group.items():
            sample = dset[:n]
            if sample.ndim == 1:
                columns[name] = sample
            elif sample.ndim == 2 and sample.shape[1] == 3:
                for i, axis in enumerate("xyz"):
                    columns[f"{name}_{axis}"] = sample[:, i]
            else:
                columns[name] = list(sample)
    return pd.DataFrame(columns)
    
    
# ---------------------------------------------------------------------------
# Notebook 2 and 3
# ---------------------------------------------------------------------------

def load_transmission(filt_path):
    """Filter CSV (wavelength_nm, R, std) -> (wave_A, resp)."""
    filt_wave_nm, filt_resp = np.loadtxt(filt_path, delimiter=',', usecols=(0, 1), unpack=True)
    return filt_wave_nm * 10, filt_resp


def get_splus_filters(base_dir):
    """Dictionary {band: path_csv} 12 S-PLUS filters."""
    return {
        'u': f'{base_dir}/uJAVA.csv',
        'J0378': f'{base_dir}/J0378.csv',
        'J0395': f'{base_dir}/J0395.csv',
        'J0410': f'{base_dir}/J0410.csv',
        'J0430': f'{base_dir}/J0430.csv',
        'g': f'{base_dir}/gSDSS.csv',
        'J0515': f'{base_dir}/J0515.csv',
        'r': f'{base_dir}/rSDSS.csv',
        'J0660': f'{base_dir}/J0660.csv',
        'i': f'{base_dir}/iSDSS.csv',
        'J0861': f'{base_dir}/J0861.csv',
        'z': f'{base_dir}/zSDSS.csv',
    }


def make_2d_header(header, subhalo_id):
    """Datacube Header -> header 2D (removes spectral axis, adds OBJECT)."""
    hdr = header.copy()
    hdr['NAXIS'] = 2
    hdr.remove('NAXIS3', ignore_missing=True)
    hdr.remove('CUNIT3', ignore_missing=True)
    hdr['OBJECT'] = subhalo_id
    return hdr


def connect_ds9():
    """Open DS9 and wait for it to register in XPA."""
    import pyds9

    subprocess.Popen(['ds9'])
    while not pyds9.ds9_targets():
        time.sleep(1)
    return pyds9.DS9()


def convert_to_jansky(data, header):
    """Convert SKIRT output (MJy/sr) to Jy/pixel using CDELT1 from header."""
    pixel_scale_arcsec = header['CDELT1']  # arcsec/pixel
    pixel_scale_rad = pixel_scale_arcsec * (np.pi / 180 / 3600)  # arcsec -> rad
    omega_pixel_sr = pixel_scale_rad**2  # sr per pixel
    return data * omega_pixel_sr * 1e6  # MJy/sr -> Jy/pixel


def convolve_filter(wave, cube, filt_path):
    filt_wave, filt_resp = load_transmission(filt_path)
    resp = np.interp(wave, filt_wave, filt_resp, left=0, right=0)
    trans = (wave * resp).astype(cube.dtype)

    mask = trans > 0
    wave_m = wave[mask]
    trans_m = trans[mask]
    cube_m = cube[mask]

    num = trapz(cube_m * trans_m[:, None, None], wave_m, axis=0)
    den = trapz(trans_m, wave_m)
    return num / den


def compute_target_pixels(z, field_phys, pixel_scale, cosmo):
    """theta = field_phys / D_A(z); num_pixels = theta / pixel_scale."""
    d_ang = cosmo.angular_diameter_distance(z)
    theta_arcsec = ((field_phys.to(u.Mpc) / d_ang) * u.rad).to(u.arcsec)
    num_pixels = int(np.round((theta_arcsec / pixel_scale).value))
    return d_ang, theta_arcsec, num_pixels


def photon_gain_factor(aperture_m, filt_path):
    """
    g = h / (A_aper * integral(R/lambda dlambda))
    eq. C.4, Zhou et al. 2025, A&A 700, A120 (arXiv:2506.15060), via Ryon 2023 ACS Handbook 9.2.
    """
    h_planck = 6.62607015e-34  # J s
    filt_wave, filt_resp = load_transmission(filt_path)
    int_R_over_lambda = trapz(filt_resp / filt_wave, filt_wave)
    area_aper_m2 = np.pi * (aperture_m / 2) ** 2
    return h_planck / (area_aper_m2 * int_R_over_lambda)


def surface_brightness_to_counts(img, pixel_scale, n_exp, t_exp, g_phys):
    """counts = N_exp * t_exp * Omega_pix * I_src / g_phys (1e-26: Jy -> SI)."""
    pixel_area_sr = pixel_scale.to(u.rad).value ** 2
    I_src_Jy = img * 1e6
    return n_exp * t_exp * pixel_area_sr * I_src_Jy * 1e-26 / g_phys


def apply_psf(image, fwhm_arcsec, pixel_scale):
    """Convolve with a gaussian PSF; sigma = FWHM / (2*sqrt(2*ln2))."""
    fwhm_pix = fwhm_arcsec / pixel_scale.to(u.arcsec).value
    sigma_pix = fwhm_pix / (2 * np.sqrt(2 * np.log(2)))
    kernel = Gaussian2DKernel(sigma_pix)
    return convolve_fft(image, kernel, normalize_kernel=True)


def apply_noise(image, n_exp, t_exp, sky_bkg, dark_current, read_noise):
    """
    Poisson(signal + sky + dark) + exposure-based gaussian readout.
    Eqs. 7-9, 11-12, Sec. 2.4, Zhou et al. 2025, A&A 700, A120 (arXiv:2506.15060).
    """
    mean_bkg = n_exp * t_exp * (sky_bkg + dark_current)
    lam = np.clip(image + mean_bkg, 0, None)
    noisy = np.random.poisson(lam).astype(np.float64)
    for _ in range(int(n_exp)):
        noisy += np.random.normal(0, read_noise, size=noisy.shape)
    return noisy - mean_bkg
