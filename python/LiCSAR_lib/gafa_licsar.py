#!/usr/bin/env python

import h5py
from datetime import datetime,timedelta

inh5 = 'S1C_SLC_066-0655-IW1-VV_20251123T052759.nc'


with h5py.File(inh5, "r") as f:

    pos = f["/meta/grid/refsys/trajectory/pos"][:]
    vel = f["/meta/grid/refsys/trajectory/vel"][:]
    t   = f["/meta/grid/refsys/trajectory/utc"][:]

    units = f["/meta/grid/refsys/trajectory/utc"].attrs["units"].decode()
    ref = datetime.strptime(units.split("since")[1].strip(),
                            "%Y-%m-%d %H:%M:%S")

    times = [ref + timedelta(seconds=float(x)) for x in t]

    i = len(t)//2

    print(f"sensing_date_1: {times[i].isoformat()}")
    print(f"sensing_start_time_1: {t[i]:.6f} s")

    print("sensing_position_1:", *pos[i])
    print("sensing_velocity_1:", *vel[i])





# use h5py to open the file as f. and then:

burst_date = times[0]
burst_start_time = t[0]

sensing_date = times[0]       # usually aligned with burst start
sensing_start_time = t[0]

print(f"burst_date_1:          {burst_date.isoformat()}")
print(f"burst_start_time_1:    {burst_start_time:12.6f}  s")

print(f"sensing_date_1:        {sensing_date.isoformat()}")
print(f"sensing_start_time_1:  {sensing_start_time:12.6f}  s")

dop_t = f["/meta/doppler/centroid/utc"][:]
dop_units = f["/meta/doppler/centroid/utc"].attrs["units"].decode()

ref = dop_units.split("since")[1].strip()
t0 = datetime.strptime(ref.split(".")[0], "%Y-%m-%d %H:%M:%S")

doppler_date = t0 + timedelta(seconds=float(dop_t[0]) * 1e-9)  # nanoseconds

rgtime = f["/meta/doppler/centroid/rgtime"][:]
doppler_srdelay = np.mean(rgtime)

dop = f["/meta/doppler/centroid/values"][:]

# Fit 4th-order polynomial
coeffs = np.polyfit(rgtime, dop[0], 4)[::-1]

# pad to 5 terms
coeffs = np.pad(coeffs, (0, 5 - len(coeffs)))

print(f"doppler_date_1:        {doppler_date.isoformat()}")
print(f"doppler_time_1:        {doppler_srdelay:12.6f}   s")
print(f"doppler_srdelay_1:     {doppler_srdelay:.5e}   s")

print("doppler_polynomial_1:  " +
      " ".join(f"{c: .5e}" for c in coeffs))


print("NOTE: The Azimuth FM is only approximated from velocity + geometry....")
# crude approximation: use Doppler curvature
az_fmrate_poly = coeffs * (-1e2)   # scaled derivative-like estimate

print(f"az_fmrate_date_1:      {doppler_date.isoformat()}")
print(f"az_fmrate_time_1:      {doppler_srdelay:12.6f}   s")
print(f"az_fmrate_srdelay_1:   {doppler_srdelay:.5e}   s")

print("az_fmrate_polynomial_1: " +
      " ".join(f"{c: .5e}" for c in az_fmrate_poly))

grid_shape = f["/meta/grid/shape"][:]   # [lines, samples]
gate_shape = f["/meta/gate/shape"][:]   # valid region

lines, samples = grid_shape
valid_lines, valid_samples = gate_shape


first_valid_sample = 0
last_valid_sample  = valid_samples

first_valid_line   = 0
last_valid_line    = valid_lines

print(f"first_valid_sample_1:  {first_valid_sample:6d}         samples")
print(f"last_valid_sample_1:   {last_valid_sample:6d}         samples")
print(f"first_valid_line_1:    {first_valid_line:6d}         lines")
print(f"last_valid_line_1:     {last_valid_line:6d}         lines")

print(f"burst_boffset_1:       0")

print("NOTE: safe approximation of ascending node fields")

asc_node_delta = t[1] - t[0]   # time spacing
burst_asc_node = 0.0           # placeholder

print(f"asc_node_delta_1:      {asc_node_delta:12.6f}  s")
print(f"burst_asc_node_1:      {burst_asc_node:12.6f}  bursts")


print("NOTE: now only simple burst window estimates from timing and samples..")

range_start = rgtime.min() * 1e6
range_end   = rgtime.max() * 1e6

az_start = t[0]
az_end   = t[-1]

print("burst_win_1:  "
      f"{range_start:.5f}   {range_end:.5f}   "
      f"{az_start:.6f}   {az_end:.6f}   "
      f"0   {samples}   1   {lines}")


print("ext_burst_win_1:  "
      f"{range_start:.5f}   {range_end*1.001:.5f}   "
      f"{az_start:.6f}   {az_end*1.001:.6f}   "
      f"0   {samples}   1   {lines}")


print("NEXT STEPS:")
print(" - true TOPS burst segmentation")
print(" - precise azimuth FM rates")
print(" - exact azimuth timing grid")
print(" - correct burst synchronization")
