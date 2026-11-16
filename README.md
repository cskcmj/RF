# Simple RF atmospheric attenuation and beam spread

`rf_atmosphere.py` calculates a Cartesian parameter sweep and writes a CSV table for UHF, L, S, C, X, K and Ka frequencies. It estimates dry-air and water-vapour absorption, cloud/fog attenuation, idealized Gaussian beam diameter, transmitted pulse energy and whole-beam peak power remaining after atmospheric attenuation.

It runs on Python 3.9 or later using the standard library. All coefficient data is bundled; calculation requires no internet connection or package installation. Keep `rf_geometry.py` and the supplied coefficient files with the script.

## Run

The repository includes the generated example table as
[`atmosphere_results.csv.gz`](https://github.com/cskcmj/RF/blob/main/atmosphere_results.csv.gz),
a lossless compressed copy of the 122,472-row CSV. Download it and extract with
`gzip -dk atmosphere_results.csv.gz`, or generate a fresh CSV with the commands
below. Its settings and model provenance are in
`atmosphere_results.csv.metadata.json`. The calculator package contains the
source code and configuration; the example table is a separate download.

From this directory:

```sh
python3 rf_atmosphere.py --output atmosphere_results.csv
python3 rf_atmosphere.py --config example_config.json --output custom_results.csv
python3 rf_atmosphere.py --geometry slant --output slant_results.csv
python3 rf_atmosphere.py --wide-max-range-mi 10 --narrow-min-range-mi 2 --output overlap_results.csv
python3 rf_atmosphere.py --config example_config.json --dry-run
python3 rf_atmosphere.py --write-example-config my_config.json
python3 rf_atmosphere.py --help
python3 -m unittest discover -v
```

`--dry-run` validates the configuration and reports the number of rows without writing output. `--write-example-config PATH` writes an editable configuration. `--geometry horizontal|slant` selects path geometry. `--wide-max-range-mi` and `--narrow-min-range-mi` override the configurable range limits for the two beam regimes, including an overlap. Every calculation CSV has a sidecar named `OUTPUT.csv.metadata.json` that records its configuration, model assumptions and provenance.

The default horizontal sweep produces **122,472 rows**: 7 frequencies × 3 altitudes × 3 density scales × 3 weather cases × 24 range/angle choices × 3 pulse durations × 3 peak powers × 3 initial diameters. There are 24 range/angle choices because the default 1-mile and 5-mile ranges each include both narrow and wide beams. This is a sample grid, not every possible value in the specified intervals.

## Beam and path interpretation

The beam has a Gaussian intensity envelope. All supplied beam opening angles are **full 1/e² divergence angles**. Initial beam size means **diameter at the source-plane Gaussian waist**, measured at the 1/e² intensity contour. The source is not focused onto a distant point. If the intended initial sizes are radii, that interpretation needs confirmation: converting them to diameters may require changing the current 0.2–2 m diameter limits.

Select **horizontal** or **slant** propagation. Horizontal is the default. Slant paths use the configured start/end altitudes and a simple vertical atmosphere profile. For a real antenna, its radiation pattern, aperture illumination, focusing and pointing information would improve the beam estimate.

Default wide beams are 15°, 22.5° and 30°, selected through 5 miles. Default narrow beams are 1°, 3° and 5°, selected from 1 mile onward. Both regimes are calculated from **1 through 5 miles, inclusive**. Change either range limit in JSON or with the CLI options above; there is no fixed short/long cutoff. The narrow minimum range must not exceed the wide maximum range, so every range has at least one regime.

## Default frequency samples

| Label | Frequency (GHz) | Accepted range (GHz) |
| --- | ---: | ---: |
| UHF | 0.5 | 0.3–1 |
| L | 1.5 | 1–2 |
| S | 3 | 2–4 |
| C | 6 | 4–8 |
| X | 10 | 8–12 |
| K | 22.23508 | 18–27 |
| Ka | 35 | 27–40 |

Labels identify configurable frequency samples, not averages over a band. The script uses microwave engineering band limits and restricts its UHF label to 0.3–1 GHz; general radio UHF terminology extends to 3 GHz and overlaps L and S. Boundary frequencies are accepted in either adjacent label. K is sampled at the 22.23508 GHz water-vapour resonance, so it is deliberately not representative of every K-band frequency. Add additional frequencies under a label to examine its variation.

## Configuration

Copy `example_config.json`, edit it and pass its path with `--config`. Units are explicit in each key. A configuration can contain only the keys to override. Each supplied key replaces the entire corresponding default value; for example, a `frequencies_ghz` object containing only `X` selects only X-band samples.

| Key | Meaning and default |
| --- | --- |
| `frequencies_ghz` | Object mapping band labels to frequency lists; defaults to the seven samples above. |
| `altitudes_ft` | Altitude above mean sea level: `[0, 5000, 10000]`. Accepted: 0–11,000 m, approximately 0–36,089.24 ft. |
| `path_geometry` | `"horizontal"` (default) or `"slant"`. In slant mode `altitudes_ft` supplies start altitudes. |
| `slant_end_altitudes_ft` | End altitudes crossed with start altitudes in slant mode: `[5000, 10000]`; same altitude bounds as the starts. |
| `slant_integration_steps` | Positive integer count of midpoint samples along a slant atmosphere profile; default `64`. |
| `density_scales` | Positive multipliers for reference moist-air mass density at each altitude, temperature and humidity: `[0.8, 1.0, 1.2]`. |
| `air_densities_kg_m3` | Optional explicit bulk moist-air density list. Default `null`; use this instead of `density_scales` when supplying measured density. |
| `temperature_k` | Optional fixed temperature for all cases, 223.15–323.15 K. Default `null` uses the ISA tropospheric temperature lapse with altitude. |
| `weather_cases` | List of weather objects, described below. |
| `ranges` | List of `{ "value": number, "unit": "ft"/"mi"/"m"/"km" }` objects. Default: 1000 ft, 3000 ft, 1 mi, 5 mi, 10 mi and 50 mi. Miles are statute miles. |
| `pulse_durations_ns` | Rectangular-envelope pulse durations: `[10, 100, 1000]` ns. Accepted: 10–1000 ns. |
| `peak_powers_w` | Launched RF peak power during each pulse: `[1000, 1000000, 1000000000]` W. Accepted: 1 kW–1 GW; this is not electrical input power or EIRP. |
| `initial_beam_diameters_m` | Source-plane 1/e² intensity diameters: `[0.2, 1.0, 2.0]` m. Accepted: 0.2–2 m. |
| `beam_quality_m2` | Minimum Gaussian beam quality factor M², a value at least 1; default `1.0`. Sets the diffraction minimum for a given frequency and diameter. |
| `divergence_mode` | `"range_presets"` (default), `"diffraction_limited"`, or `"specified"`; described below. |
| `wide_beam_max_range_mi` | Wide beams are selected at or below this range; default `5.0` statute miles. |
| `narrow_beam_min_range_mi` | Narrow beams are selected at or above this range; default `1.0` statute miles. Must not exceed `wide_beam_max_range_mi`. |
| `narrow_full_angles_deg` | Full opening angles for the narrow regime: `[1, 3, 5]` degrees. |
| `wide_full_angles_deg` | Full opening angles for the wide regime: `[15, 22.5, 30]` degrees. |
| `specified_full_angles_deg` | Full opening angles applied at every range in `"specified"` mode: `[3]` degrees. |
| `max_paraxial_half_angle_rad` | Heuristic validity threshold for the Gaussian model; default `0.3` rad. Accepted: greater than 0, at most 0.3 rad. |

`density_scales` and `air_densities_kg_m3` are alternative ways to set bulk density. When using an explicit density list, set `density_scales` to `null`. All list dimensions are crossed with one another, so table size grows rapidly when adding samples.

Each weather case has four fields:

```json
{
  "name": "medium_fog",
  "relative_humidity_pct": 100,
  "liquid_water_g_m3": 0.05,
  "cloud_path_fraction": 1.0
}
```

Defaults are `clear` (50% relative humidity, no liquid water), `medium_fog` (100%, 0.05 g/m³), and `thick_fog` (100%, 0.5 g/m³). All have `cloud_path_fraction: 1.0`. The liquid-water values correspond to illustrative medium and thick fog in ITU-R P.840; they do not describe all clouds or fogs. `cloud_path_fraction` is the fraction of the propagation distance inside the specified liquid cloud or fog. Gas absorption applies along the entire path. A fraction of 0.2 therefore applies cloud/fog attenuation along 20% of the range.

Bulk air mass density (kg/m³), water-vapour density (g/m³) and suspended liquid-water content (g/m³) are separate quantities. Humidity and temperature determine vapour pressure; the program derives consistent dry pressure from the specified bulk density. The output includes these atmospheric quantities so each case can be inspected.

## Path geometry

In horizontal mode, each range lies at one constant altitude and `end_altitude_ft` equals `altitude_ft`. In slant mode, start altitudes, end altitudes and ranges are crossed. **Range is the straight path length**, not horizontal ground distance. Elevation is `asin((end_altitude_m - altitude_m) / range_m)`, positive for ascending paths and negative for descending paths. A zero-length path is permitted only between equal altitudes.

Combinations with a path shorter than the endpoint altitude difference are skipped. The metadata's `rows_skipped_impossible_geometry` counts skipped **full output rows**, including frequency, weather, beam and pulse combinations. It does not count only unique endpoint/range combinations.

Along a slant path, altitude varies linearly with path distance. Midpoint integration evaluates temperature, humidity-dependent vapour density, air density and gas/cloud loss at each sample altitude. With density scales, the scale multiplies local reference mass density throughout the profile. Explicit `air_densities_kg_m3` instead imposes the supplied mass density throughout the path, while local pressure is recalculated consistently. The optional fixed temperature also applies throughout the path.

Relative humidity and liquid-water content stay at their weather-case values along the path. `cloud_path_fraction` applies uniform fractional cloud occupancy along this profile; it does not define a cloud layer at a particular altitude. The CSV's specific-attenuation columns contain path-mean rates. Source atmospheric properties remain in the existing temperature/density/pressure columns, while `path_mean_temperature_k` and `path_mean_air_density_kg_m3` report the path averages.

The default slant sweep produces **193,914 rows** and skips **51,030 impossible-geometry rows**. Same-height slant cases are retained and agree with horizontal calculations.

## Models and equations

### Atmosphere and absorption

Altitude supplies an ISA tropospheric temperature and pressure reference, treating geometric altitude as the standard altitude approximation over 0–11 km. The optional temperature and density inputs allow departures from that reference. The updated Buck saturation-vapour-pressure expression over liquid water, without its pressure enhancement factor, and relative humidity determine vapour partial pressure e. Relative humidity is therefore defined with respect to liquid water even below freezing; liquid droplets at those temperatures are treated as supercooled, and the output is flagged.

Reference moist-air mass density uses the ISA total pressure at the chosen altitude, the selected temperature and relative humidity. `density_scales` multiplies that reference mass density. With target mass density and vapour density in kg/m³:

```text
vapour_density_kg_m3 = e_hpa * 100 / (R_vapour * temperature_k)
dry_pressure_hpa = (target_air_density_kg_m3 - vapour_density_kg_m3) * R_dry * temperature_k / 100
total_pressure_hpa = dry_pressure_hpa + e_hpa
```

For explicit `air_densities_kg_m3`, the reported `density_scale` is the supplied mass density divided by that reference density. Configurations producing a nonpositive dry pressure are rejected. For the attenuation equations, vapour density is calculated using the ITU conversion ρᵥ = 216.7e/T g/m³, with pressure in hPa and temperature in K. The gas model uses **dry pressure**, p = P_total − e, rather than passing total pressure as dry pressure.

Gas attenuation follows the line-by-line equations in ITU-R P.676-13, Annex 1, including the dry continuum and the water-vapour pseudo-line. The local 44-line oxygen and 35-line water tables agree numerically with the recommendation's Tables 1 and 2. Cloud/fog attenuation uses ITU-R P.840-9's double-Debye liquid-water coefficient (K_l(f,T)).

For range r in km and cloud-path fraction q, use the local rate for a horizontal path or its length-weighted mean for a slant path:

```text
A_dry_dB   = gamma_dry_dB_per_km * r
A_vapour_dB = gamma_vapour_dB_per_km * r
A_cloud_dB = Kl_dB_per_km_per_g_m3 * liquid_water_g_m3 * q * r
A_total_dB = A_dry_dB + A_vapour_dB + A_cloud_dB
transmission = 10 ** (-A_total_dB / 10)
```

The recommendation's stated gas-model scope begins at 1 GHz. The same analytical expression is evaluated below 1 GHz to include UHF, and those rows are explicitly flagged as extrapolations.

### Beam diameter and spreading

`divergence_mode` selects how the Gaussian-equivalent divergence is chosen:

- `range_presets` uses the wide and narrow angle lists with their overlapping range limits. Overlap rows retain distinct `wide` and `narrow` regime labels.
- `specified` applies every angle in `specified_full_angles_deg` at every range.
- `diffraction_limited` uses the minimum divergence for the configured M², diameter and wavelength, with one angle choice per range.

Let D₀ be the initial diameter, w₀ = D₀/2, B = configured minimum M², and λ = c/f:

```text
minimum_half_angle_rad = B * wavelength_m / (pi * w0_m)
used_half_angle_rad = radians(requested_full_divergence_deg / 2)
effective_beam_quality_m2 = used_half_angle_rad * pi * w0_m / wavelength_m
rayleigh_range_m = w0_m / used_half_angle_rad
w_m(r) = w0_m * sqrt(1 + (range_m / rayleigh_range_m)**2)
final_diameter_m = 2 * w_m(r)
area_ratio = (w_m(r) / w0_m)**2
on_axis_spreading_loss_dB = 10 * log10(area_ratio)
```

In `diffraction_limited` mode, `used_half_angle_rad` is the minimum half-angle above. In the other modes, it is the requested half-angle. These use the Gaussian paraxial angle convention: angle in radians approximates the far-field radius-growth slope. The difference between this and a trigonometric opening angle is about 2.35% in far-field width at a full 30° angle.

The spreading quantity describes the change in ideal Gaussian on-axis intensity caused by growth in beam area. It does not remove energy from the entire beam. Atmospheric attenuation and area dilution are reported separately; free-space path loss must not be applied again to the same whole-beam power calculation.

A requested angle below the diffraction minimum receives `below_diffraction_limit`; its final diameter and dependent beam-intensity quantities are blank. This occurs frequently with narrow UHF/L beams at the requested diameters. For example, at 0.5 GHz a 2 m diameter and M² = 1 already imply a minimum full angle of about 21.87°, so a 1–5° request is not feasible in this model. The supplied angle is never silently widened.

The Gaussian formula also assumes paraxial propagation. A used half-angle above `max_paraxial_half_angle_rad` receives `outside_paraxial_limit`, with beam predictions blank. The 0.3 rad threshold is a documented screening choice, not a guarantee of antenna-model accuracy. The default 15–30° full wide-beam angles have 7.5–15° half-angles and pass that screen, while remaining approximate. `--allow-invalid-beam-estimates` includes diagnostic numerical beam values for invalid requests and prefixes the status with `unreliable_estimate_`; it does not make a below-diffraction request physically feasible.

When the requested angle and initial diameter are fixed, the calculated final diameter is independent of frequency for feasible rows. Frequency still changes atmospheric attenuation, the diffraction minimum and the required effective M².

### Pulses and power

```text
input_pulse_energy_j = peak_power_w * pulse_duration_s
transmitted_peak_power_w = peak_power_w * transmission
transmitted_pulse_energy_j = input_pulse_energy_j * transmission
on_axis_peak_power_density_w_m2 = 2 * transmitted_peak_power_w / (pi * w_m(r)**2)
pulse_length_m = c * pulse_duration_s
first_null_offset_mhz = 1 / pulse_duration_s / 1e6
first_null_offset_fraction_of_carrier = 1 / (pulse_duration_s * frequency_hz)
```

The power after attenuation is integrated across the whole propagated beam. It is **not received power**: receiver aperture, pointing, efficiency and antenna pattern have not been supplied. A pulse repetition rate is needed to calculate average power; none is assumed.

Pulse duration changes pulse energy and the characteristic bandwidth. Peak power scales power, energy and intensity. Neither changes attenuation or beam diameter in this linear model. For an ideal rectangular pulse, 1/τ is the spectral offset from the carrier to its first null; the null-to-null main-lobe width is 2/τ, and sidelobes extend farther. It is not a hard occupied-bandwidth cutoff. A first-null offset greater than 10% of the carrier produces `large_pulse_spectral_offset` in `model_flags`, because the model evaluates loss and diffraction at the carrier only. This is a screening heuristic. No pulse dispersion or temporal broadening is calculated. `pulse_length_m` uses the vacuum speed of light.

## Reading the output

The CSV has these 56 columns. `final_beam_area_m2` is the area inside the Gaussian 1/e² diameter; the Gaussian distribution extends beyond that contour. `cloud_db_per_km` is the path-mean in-cloud rate, before applying `cloud_path_fraction`.

| Columns | Meaning |
| --- | --- |
| `band`, `frequency_ghz`, `wavelength_m` | Frequency sample and free-space wavelength. |
| `range_m`, `range_ft`, `range_mi` | Same one-way propagation distance in three units. |
| `path_geometry`, `altitude_ft`, `altitude_m`, `end_altitude_ft`, `end_altitude_m`, `elevation_deg` | Geometry, start/end altitudes and signed elevation angle. |
| `temperature_k`, `path_mean_temperature_k`, `path_mean_air_density_kg_m3` | Source temperature and path-mean atmospheric conditions. |
| `density_scale`, `air_density_kg_m3`, `total_pressure_hpa`, `dry_pressure_hpa` | Air state derived consistently from the configured density. |
| `weather`, `relative_humidity_pct`, `water_vapour_density_g_m3` | Weather label, liquid-water-relative humidity and vapour mass density. |
| `liquid_water_g_m3`, `cloud_path_fraction`, `cloud_path_km` | Liquid-droplet content and length of cloud/fog along the path. |
| `dry_gas_db_per_km`, `water_vapour_db_per_km`, `cloud_db_per_km` | Path-mean specific attenuation contributions. |
| `gas_attenuation_db`, `cloud_attenuation_db`, `atmospheric_attenuation_db`, `atmospheric_power_transmission` | Integrated path losses and remaining power fraction. |
| `initial_beam_diameter_m`, `beam_quality_m2`, `rayleigh_range_m`, `nominal_divergence_half_angle_rad`, `beam_model_status` | Source waist, nominal Gaussian parameters and validity status. |
| `divergence_regime`, `requested_full_divergence_deg`, `used_full_divergence_deg`, `diffraction_min_full_divergence_deg`, `effective_beam_quality_m2` | Regime label, full requested/used/minimum angles and implied beam quality. |
| `final_beam_diameter_m`, `final_beam_area_m2`, `beam_area_increase_factor`, `on_axis_spreading_loss_db` | Beam size and geometric on-axis intensity dilution. |
| `pulse_duration_ns`, `pulse_length_m`, `first_null_offset_mhz`, `first_null_offset_fraction_of_carrier` | Pulse duration, spatial length and spectral screening quantities. |
| `peak_power_w`, `input_pulse_energy_j`, `transmitted_peak_power_w`, `transmitted_pulse_energy_j`, `on_axis_peak_power_density_w_m2` | Input and surviving whole-beam power/energy, plus Gaussian on-axis peak intensity. |
| `model_flags` | Semicolon-separated assumptions and validity flags. |

`beam_model_status` is `paraxial_estimate` for a feasible row within the chosen angle limit, `below_diffraction_limit` or `outside_paraxial_limit` when beam predictions are blank, or the corresponding status prefixed with `unreliable_estimate_` when invalid estimates were requested explicitly. Blank beam values are not zero-size beams. In `diffraction_limited` mode, `requested_full_divergence_deg` is blank because no angle was requested; `used_full_divergence_deg` reports the calculated minimum. `nominal_divergence_half_angle_rad` is the half-angle used in the calculation. Inspect `beam_model_status` and `model_flags` before comparing numerical results.

Every row includes `linear_model_only` and `carrier_frequency_only`, plus either `horizontal_uniform_atmosphere` or `slant_profile_midpoint_integration`. Conditional flags are `itu_models_extrapolated_below_1ghz`, `large_pulse_spectral_offset`, `supercooled_liquid_assumed_no_ice`, `requested_divergence_below_diffraction_limit` and `beam_outside_paraxial_limit`.

The metadata sidecar records the effective configuration, model version, model/source references, coefficient SHA-256 hashes, output row count, beam-validity counts and skipped impossible-geometry rows. Retain it with the CSV when sharing results.

## Limits

This is a first-order, one-way calculation on a straight horizontal or slant path. Horizontal conditions are uniform; slant profiles use simple altitude-dependent reference conditions. It includes no terrain, obstruction, Earth curvature, horizon check, ground reflection, multipath, rain, snow, ice-cloud attenuation, turbulence, ducting, refraction, antenna sidelobes or pointing error. A long near-ground path is not necessarily unobstructed simply because a CSV row exists.

Liquid clouds and fog are represented by their liquid-water content, not visibility or droplet-resolved scattering. A uniform attenuating medium changes beam amplitude, not the normalized transverse beam shape; atmosphere-dependent beam widening would require an additional refractive/turbulence model.

The peak-power sweep extends to 1 GW as requested, but these calculations assume linear propagation. They do not predict electrical breakdown, ionization, plasma formation, heating or other nonlinear changes to the atmosphere. A numerical row at high power does not establish that those assumptions hold.

## References and provenance

- [ITU-R P.676-13: gas absorption, Annex 1](https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.676-13-202208-I%21%21PDF-E.pdf).
- [ITU-R P.840-9: attenuation due to clouds and fog](https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.840-9-202308-I%21%21PDF-E.pdf).
- [ITU-R P.835: reference atmospheres](https://www.itu.int/rec/R-REC-P.835/en).
- [Buck Research CR-1A operating manual, Appendix 1: humidity conversion equations](https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf).
- [ITU-Rpy oxygen coefficient transcription](https://raw.githubusercontent.com/inigodelportillo/ITU-Rpy/master/itur/data/676/v12_lines_oxygen.txt) and [water-vapour coefficient transcription](https://raw.githubusercontent.com/inigodelportillo/ITU-Rpy/master/itur/data/676/v12_lines_water_vapour.txt), checked against P.676-13 Tables 1 and 2.
- [ITU-Rpy gas-model validation examples](https://itu-rpy.readthedocs.io/en/latest/validation/itur676_12.html), used for independent numerical reference checks.
- [Edmund Optics: Gaussian beam propagation](https://www.edmundoptics.ca/knowledge-center/application-notes/lasers/gaussian-beam-propagation/).
