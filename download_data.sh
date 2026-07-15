#!/usr/bin/env bash
# Download the TNG50 cutouts and SKIRT datacubes from the GitHub Release.
#
# These files are attached to the Release, not committed to the repo —
# the (gzip-compressed) datacubes are up to ~600 MB and the cutouts ~150 MB,
# too much for git.

set -e

# ---------------------------------------------------------------------
# Your repo and the release tag
# ---------------------------------------------------------------------
REPO="GoncalvesGus/SPLUS_SKIRT_tutorial"
TAG="v1.0-data"
# ---------------------------------------------------------------------

BASE="https://github.com/${REPO}/releases/download/${TAG}"

download() {
    local dir="$1"
    shift
    mkdir -p "$dir"
    cd "$dir"
    for f in "$@"; do
        if [ -f "$f" ]; then
            echo "  [have]  $f"
        else
            echo "  [get ]  $f"
            curl -fL --progress-bar -O "${BASE}/${f}"
        fi
    done
    cd ..
}


SKIRT_FILES=(
    "10_splus_cube_total.fits.gz"
    "406122_splus_cube_total.fits.gz"
)

TNG_FILES=(
    "cutout_TNG50-1_snap97_subhalo10.hdf5"
    "subhalo10_firefly_positions.hdf5"
    "cutout_TNG50-1_snap98_subhalo406122.hdf5"
    "subhalo406122_firefly_positions.hdf5"
)

download "SKIRT_datacubes" "${SKIRT_FILES[@]}"
download "TNG50_hdf5_cutouts" "${TNG_FILES[@]}"

echo
echo "Done. Files are in place under SKIRT_datacubes/ and TNG50_hdf5_cutouts/."
