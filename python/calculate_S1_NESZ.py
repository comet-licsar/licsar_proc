#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import numpy as np
import sys

if len(sys.argv) != 3:
    print("Usage:")
    print("  python calc_nesz.py noise.xml calibration.xml")
    sys.exit(1)

noise_file = sys.argv[1]
cal_file = sys.argv[2]

# ---------------------------------------------------------------------
# Read noise file
# ---------------------------------------------------------------------
noise_root = ET.parse(noise_file).getroot()

nrv = noise_root.find(".//noiseRangeVector")
nav = noise_root.find(".//noiseAzimuthVector")

pixels_noise = np.array(
    list(map(float, nrv.find("pixel").text.split()))
)

noise_range = np.array(
    list(map(float, nrv.find("noiseRangeLut").text.split()))
)

noise_azimuth = np.array(
    list(map(float, nav.find("noiseAzimuthLut").text.split()))
)

# Usually only one azimuth profile exists
# Use mean azimuth correction for a quick estimate
noise_az = np.mean(noise_azimuth)

# ---------------------------------------------------------------------
# Read calibration file
# ---------------------------------------------------------------------
cal_root = ET.parse(cal_file).getroot()

cv = cal_root.find(".//calibrationVector")

pixels_cal = np.array(
    list(map(float, cv.find("pixel").text.split()))
)

sigma_nought = np.array(
    list(map(float, cv.find("sigmaNought").text.split()))
)

# ---------------------------------------------------------------------
# Ensure same grid
# ---------------------------------------------------------------------
if not np.array_equal(pixels_noise, pixels_cal):
    sigma_nought = np.interp(
        pixels_noise,
        pixels_cal,
        sigma_nought
    )

# ---------------------------------------------------------------------
# NESZ calculation
# ---------------------------------------------------------------------
noise_power = noise_range * noise_az

nesz_linear = noise_power / sigma_nought**2

nesz_db = 10.0 * np.log10(nesz_linear)

# ---------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------
print("")
print(f"Mean azimuth noise factor : {noise_az:.6f}")
print(f"NESZ min                 : {nesz_db.min():.2f} dB")
print(f"NESZ max                 : {nesz_db.max():.2f} dB")
print(f"NESZ mean                : {nesz_db.mean():.2f} dB")
print("")

print("Pixel     NESZ[dB]")
for p, n in zip(pixels_noise[::20], nesz_db[::20]):
    print(f"{int(p):6d}   {n:8.2f}")