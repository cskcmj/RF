#!/usr/bin/env python3
"""Sweep simple horizontal/slant RF propagation cases into CSV (Python 3.9+, no deps).

Gas: ITU-R P.676-13 Annex 1. Cloud liquid water: P.840-9 Annex 1, sections 1–2.
Beam: source-plane Gaussian waist, 1/e² intensity diameter, paraxial diffraction.
See README.md for assumptions, sources, configuration, and column definitions.
"""

import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
from functools import lru_cache
from pathlib import Path
import sys
import tempfile

from rf_geometry import elevation_deg, integrate_path, valid_path

C = 299_792_458.0
FT = 0.3048
MI = 1609.344
R_DRY = 287.05
R_VAPOUR = 461.5
DATA_DIR = Path(__file__).resolve().parent / "data"
MODEL_VERSION = "1.1"

# Microwave engineering convention. UHF is restricted to 0.3–1 GHz here to
# avoid overlapping L and S; general radio UHF terminology can extend to 3 GHz.
BAND_LIMITS_GHZ = {
    "UHF": (0.3, 1.0), "L": (1.0, 2.0), "S": (2.0, 4.0),
    "C": (4.0, 8.0), "X": (8.0, 12.0), "K": (18.0, 27.0),
    "Ka": (27.0, 40.0),
}
DEFAULT_CONFIG = {
    "frequencies_ghz": {
        "UHF": [0.5], "L": [1.5], "S": [3.0], "C": [6.0],
        "X": [10.0], "K": [22.23508], "Ka": [35.0],
    },
    "altitudes_ft": [0, 5000, 10000],
    "path_geometry": "horizontal",
    "slant_end_altitudes_ft": [5000, 10000],
    "slant_integration_steps": 64,
    "density_scales": [0.8, 1.0, 1.2],
    "air_densities_kg_m3": None,
    "temperature_k": None,
    "weather_cases": [
        {"name": "clear", "relative_humidity_pct": 50,
         "liquid_water_g_m3": 0, "cloud_path_fraction": 1},
        {"name": "medium_fog", "relative_humidity_pct": 100,
         "liquid_water_g_m3": 0.05, "cloud_path_fraction": 1},
        {"name": "thick_fog", "relative_humidity_pct": 100,
         "liquid_water_g_m3": 0.5, "cloud_path_fraction": 1},
    ],
    "ranges": [{"value": 1000, "unit": "ft"}, {"value": 3000, "unit": "ft"},
               {"value": 1, "unit": "mi"}, {"value": 5, "unit": "mi"},
               {"value": 10, "unit": "mi"}, {"value": 50, "unit": "mi"}],
    "pulse_durations_ns": [10, 100, 1000],
    "peak_powers_w": [1e3, 1e6, 1e9],
    "initial_beam_diameters_m": [0.2, 1.0, 2.0],
    "divergence_mode": "range_presets",
    "wide_beam_max_range_mi": 5.0,
    "narrow_beam_min_range_mi": 1.0,
    "narrow_full_angles_deg": [1.0, 3.0, 5.0],
    "wide_full_angles_deg": [15.0, 22.5, 30.0],
    "specified_full_angles_deg": [3.0],
    "beam_quality_m2": 1.0,
    "max_paraxial_half_angle_rad": 0.3,
}

FIELDS = [
    "band", "frequency_ghz", "wavelength_m", "range_m", "range_ft", "range_mi",
    "altitude_ft", "altitude_m", "path_geometry", "end_altitude_ft", "end_altitude_m",
    "elevation_deg", "temperature_k", "path_mean_temperature_k",
    "path_mean_air_density_kg_m3", "density_scale",
    "air_density_kg_m3", "total_pressure_hpa", "dry_pressure_hpa",
    "weather", "relative_humidity_pct", "water_vapour_density_g_m3",
    "liquid_water_g_m3", "cloud_path_fraction", "cloud_path_km",
    "dry_gas_db_per_km", "water_vapour_db_per_km", "cloud_db_per_km",
    "gas_attenuation_db", "cloud_attenuation_db", "atmospheric_attenuation_db",
    "atmospheric_power_transmission", "initial_beam_diameter_m", "beam_quality_m2",
    "divergence_regime", "requested_full_divergence_deg", "used_full_divergence_deg",
    "diffraction_min_full_divergence_deg", "effective_beam_quality_m2",
    "rayleigh_range_m", "nominal_divergence_half_angle_rad", "beam_model_status",
    "final_beam_diameter_m", "final_beam_area_m2", "beam_area_increase_factor",
    "on_axis_spreading_loss_db", "pulse_duration_ns", "pulse_length_m",
    "first_null_offset_mhz", "first_null_offset_fraction_of_carrier",
    "peak_power_w", "input_pulse_energy_j", "transmitted_peak_power_w",
    "transmitted_pulse_energy_j", "on_axis_peak_power_density_w_m2",
    "model_flags",
]


def number(value, name, lower=None, upper=None, strict_lower=False):
    """Validate JSON numeric inputs; reject booleans, NaN and infinity."""
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError(f"{name} must be a number")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if lower is not None and (value < lower or (strict_lower and value == lower)):
        raise ValueError(f"{name} must be {'>' if strict_lower else '>='} {lower}")
    if upper is not None and value > upper:
        raise ValueError(f"{name} must be <= {upper}")
    return float(value)


def numeric_list(values, name, lower=0, upper=None, strict_lower=True):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} must be a nonempty list")
    return [number(v, name, lower, upper, strict_lower) for v in values]


def validate_config(config):
    unknown = set(config) - set(DEFAULT_CONFIG)
    if unknown:
        raise ValueError(f"Unknown configuration keys: {', '.join(sorted(unknown))}")
    frequencies = config["frequencies_ghz"]
    if not isinstance(frequencies, dict) or not frequencies:
        raise ValueError("frequencies_ghz must be a nonempty band-to-list mapping")
    for band, values in frequencies.items():
        if band not in BAND_LIMITS_GHZ:
            raise ValueError(f"Unknown band: {band}")
        lo, hi = BAND_LIMITS_GHZ[band]
        numeric_list(values, f"frequencies_ghz.{band}", lo, hi, False)
    # Simple troposphere only, using geometric altitude as an approximation.
    numeric_list(config["altitudes_ft"], "altitudes_ft", 0, 11000 / FT, False)
    if config["path_geometry"] not in ("horizontal", "slant"):
        raise ValueError("path_geometry must be horizontal or slant")
    numeric_list(config["slant_end_altitudes_ft"], "slant_end_altitudes_ft", 0, 11000 / FT, False)
    steps = config["slant_integration_steps"]
    if isinstance(steps, bool) or not isinstance(steps, int) or not 1 <= steps <= 10000:
        raise ValueError("slant_integration_steps must be an integer from 1 to 10000")
    if config["air_densities_kg_m3"] is not None:
        numeric_list(config["air_densities_kg_m3"], "air_densities_kg_m3")
        if config["density_scales"] is not None:
            raise ValueError("Set density_scales to null when using air_densities_kg_m3")
    else:
        numeric_list(config["density_scales"], "density_scales")
    if config["temperature_k"] is not None:
        number(config["temperature_k"], "temperature_k", 223.15, 323.15)
    numeric_list(config["pulse_durations_ns"], "pulse_durations_ns", 10, 1000, False)
    numeric_list(config["peak_powers_w"], "peak_powers_w", 1e3, 1e9, False)
    numeric_list(config["initial_beam_diameters_m"], "initial_beam_diameters_m", 0.2, 2, False)
    number(config["beam_quality_m2"], "beam_quality_m2", 1)
    number(config["max_paraxial_half_angle_rad"], "max_paraxial_half_angle_rad", 0, 0.3, True)
    if config["divergence_mode"] not in ("range_presets", "diffraction_limited", "specified"):
        raise ValueError("divergence_mode must be range_presets, diffraction_limited or specified")
    numeric_list(config["narrow_full_angles_deg"], "narrow_full_angles_deg", 1, 5, False)
    numeric_list(config["wide_full_angles_deg"], "wide_full_angles_deg", 15, 30, False)
    numeric_list(config["specified_full_angles_deg"], "specified_full_angles_deg", 0, 180, True)
    number(config["wide_beam_max_range_mi"], "wide_beam_max_range_mi", 0)
    number(config["narrow_beam_min_range_mi"], "narrow_beam_min_range_mi", 0)
    if config["narrow_beam_min_range_mi"] > config["wide_beam_max_range_mi"]:
        raise ValueError("narrow_beam_min_range_mi must be <= wide_beam_max_range_mi; overlaps are allowed")
    if not isinstance(config["ranges"], list) or not config["ranges"]:
        raise ValueError("ranges must be a nonempty list")
    for r in config["ranges"]:
        if not isinstance(r, dict) or set(r) != {"value", "unit"}:
            raise ValueError("Each range must contain exactly value and unit")
        range_metres(r)
    weather_cases = config["weather_cases"]
    if not isinstance(weather_cases, list) or not weather_cases:
        raise ValueError("weather_cases must be a nonempty list")
    names = set()
    expected = {"name", "relative_humidity_pct", "liquid_water_g_m3", "cloud_path_fraction"}
    for weather in weather_cases:
        if not isinstance(weather, dict) or set(weather) != expected:
            raise ValueError(f"Each weather case must contain exactly {sorted(expected)}")
        name = weather["name"]
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Weather names must be unique nonempty strings")
        if name[0] in "=+-@" or any(ord(c) < 32 for c in name):
            raise ValueError("Weather names cannot start with spreadsheet formula characters or contain controls")
        names.add(name)
        number(weather["relative_humidity_pct"], "relative_humidity_pct", 0, 100)
        number(weather["liquid_water_g_m3"], "liquid_water_g_m3", 0)
        number(weather["cloud_path_fraction"], "cloud_path_fraction", 0, 1)
    # Validate every thermodynamic combination before creating an output file.
    densities = config["air_densities_kg_m3"] or config["density_scales"]
    endpoints = config["altitudes_ft"] + (config["slant_end_altitudes_ft"] if config["path_geometry"] == "slant" else [])
    for altitude, density, weather in itertools.product(endpoints, densities, weather_cases):
        atmosphere(altitude, density, weather, config)
    if row_count(config) == 0:
        raise ValueError("No feasible paths: slant range must be at least the endpoint altitude difference")


def range_metres(spec):
    factors = {"m": 1.0, "km": 1000.0, "ft": FT, "mi": MI}
    if spec["unit"] not in factors:
        raise ValueError("Range unit must be m, km, ft, or mi (statute miles)")
    distance = number(spec["value"], "range value", 0) * factors[spec["unit"]]
    if not math.isfinite(distance):
        raise ValueError("Range overflows after conversion to metres")
    return distance


def standard_atmosphere(altitude_m):
    """ISA tropospheric lapse rate; geometric altitude approximation, 0–11 km."""
    number(altitude_m, "altitude_m", 0, 11000)
    temperature = 288.15 - 0.0065 * altitude_m
    pressure = 1013.25 * (temperature / 288.15) ** (9.80665 / (R_DRY * 0.0065))
    return temperature, pressure


def atmosphere(altitude_ft, density, weather, config):
    temperature_isa, pressure_isa = standard_atmosphere(altitude_ft * FT)
    temperature = config["temperature_k"] if config["temperature_k"] is not None else temperature_isa
    temp_c = temperature - 273.15
    # Updated Buck (1996) saturation pressure over liquid water, hPa.
    # The small pressure enhancement factor is omitted in this simple model.
    saturation = 6.1121 * math.exp((18.678 - temp_c / 234.5) * temp_c / (257.14 + temp_c))
    e = saturation * weather["relative_humidity_pct"] / 100.0
    vapour_density_kg = e * 100 / (R_VAPOUR * temperature)
    reference_density = (pressure_isa - e) * 100 / (R_DRY * temperature) + vapour_density_kg
    if reference_density <= 0 or pressure_isa <= e:
        raise ValueError("Temperature/humidity combination exceeds the reference total pressure")
    if config["air_densities_kg_m3"] is None:
        target_density = reference_density * density
        scale = density
    else:
        target_density = density
        scale = density / reference_density
    dry_pressure = (target_density - vapour_density_kg) * R_DRY * temperature / 100
    total_pressure = dry_pressure + e
    if dry_pressure <= 0:
        raise ValueError("Density/temperature/humidity combination gives nonpositive dry-air pressure")
    air_density = dry_pressure * 100 / (R_DRY * temperature) + vapour_density_kg
    # Use the ITU conversion exactly for consistency with its spectral equations.
    vapour_g = 216.7 * e / temperature
    return {
        "altitude_ft": altitude_ft, "altitude_m": altitude_ft * FT,
        "temperature_k": temperature, "density_scale": scale,
        "air_density_kg_m3": air_density, "total_pressure_hpa": total_pressure,
        "dry_pressure_hpa": dry_pressure, "weather": weather["name"],
        "relative_humidity_pct": weather["relative_humidity_pct"],
        "water_vapour_density_g_m3": vapour_g,
        "liquid_water_g_m3": weather["liquid_water_g_m3"],
        "cloud_path_fraction": weather["cloud_path_fraction"],
    }


@lru_cache(maxsize=2)
def spectral_lines(filename):
    with (DATA_DIR / filename).open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader)
        rows = tuple(tuple(float(v) for v in row) for row in reader if row)
    expected_count = 44 if filename == "oxygen_lines.csv" else 35
    if len(rows) != expected_count or any(len(row) != 7 for row in rows):
        raise ValueError(f"Invalid spectroscopic coefficient file: {filename}")
    return rows


def gas_specific_attenuation(f_ghz, dry_pressure_hpa, vapour_g_m3, temperature_k):
    """Return dry and wet specific loss (dB/km), P.676-13 Annex 1 eqs. 1–9.

    Input pressure is DRY pressure. Below 1 GHz is an explicit extrapolation.
    """
    f, p, rho, t = f_ghz, dry_pressure_hpa, vapour_g_m3, temperature_k
    theta = 300 / t
    e = rho * t / 216.7
    oxygen = []
    for fi, a1, a2, a3, a4, a5, a6 in spectral_lines("oxygen_lines.csv"):
        strength = a1 * 1e-7 * p * theta**3 * math.exp(a2 * (1 - theta))
        width = math.sqrt((a3 * 1e-4 * (p * theta**(0.8 - a4) + 1.1 * e * theta))**2 + 2.25e-6)
        delta = (a5 + a6 * theta) * 1e-4 * (p + e) * theta**0.8
        shape = f / fi * ((width - delta * (fi - f)) / ((fi - f)**2 + width**2)
                          + (width - delta * (fi + f)) / ((fi + f)**2 + width**2))
        oxygen.append(strength * shape)
    d = 5.6e-4 * (p + e) * theta**0.8
    continuum = f * p * theta**2 * (
        6.14e-5 / (d * (1 + (f / d)**2))
        + 1.4e-12 * p * theta**1.5 / (1 + 1.9e-5 * f**1.5))
    water = []
    for fi, b1, b2, b3, b4, b5, b6 in spectral_lines("water_vapour_lines.csv"):
        strength = b1 * 1e-1 * e * theta**3.5 * math.exp(b2 * (1 - theta))
        width = b3 * 1e-4 * (p * theta**b4 + b5 * e * theta**b6)
        width = 0.535 * width + math.sqrt(0.217 * width**2 + 2.1316e-12 * fi**2 / theta)
        shape = f / fi * (width / ((fi - f)**2 + width**2) + width / ((fi + f)**2 + width**2))
        water.append(strength * shape)
    return 0.1820 * f * (math.fsum(oxygen) + continuum), 0.1820 * f * math.fsum(water)


def cloud_specific_attenuation(f_ghz, temperature_k, liquid_water_g_m3):
    """Local liquid-droplet loss (dB/km), P.840-9 eqs. 1–10."""
    theta = 300 / temperature_k
    epsilon0 = 77.66 + 103.3 * (theta - 1)
    epsilon1, epsilon2 = 0.0671 * epsilon0, 3.52
    fp = 20.20 - 146 * (theta - 1) + 316 * (theta - 1)**2
    fs = 39.8 * fp
    real = ((epsilon0 - epsilon1) / (1 + (f_ghz / fp)**2)
            + (epsilon1 - epsilon2) / (1 + (f_ghz / fs)**2) + epsilon2)
    imag = (f_ghz * (epsilon0 - epsilon1) / (fp * (1 + (f_ghz / fp)**2))
            + f_ghz * (epsilon1 - epsilon2) / (fs * (1 + (f_ghz / fs)**2)))
    eta = (2 + real) / imag
    return 0.819 * f_ghz / (imag * (1 + eta**2)) * liquid_water_g_m3


def gaussian_beam(f_ghz, initial_diameter_m, distance_m, m2=1.0, max_angle=0.3,
                  allow_invalid=False, full_angle_deg=None):
    """Gaussian waist propagation, with an optional full divergence angle.

    The selected angle must meet the diffraction floor for the supplied waist.
    Angles in radians are the paraxial slope convention, including at 30° full.
    """
    wavelength = C / (f_ghz * 1e9)
    w0 = initial_diameter_m / 2
    minimum_theta = m2 * wavelength / (math.pi * w0)
    theta = minimum_theta if full_angle_deg is None else math.radians(full_angle_deg) / 2
    z_rayleigh = w0 / theta
    below_floor = theta < minimum_theta * (1 - 1e-12)
    status = ("below_diffraction_limit" if below_floor else
              "outside_paraxial_limit" if theta > max_angle else "paraxial_estimate")
    valid = status == "paraxial_estimate"
    result = {
        "initial_beam_diameter_m": initial_diameter_m, "beam_quality_m2": m2,
        "requested_full_divergence_deg": full_angle_deg,
        "used_full_divergence_deg": math.degrees(2 * theta),
        "diffraction_min_full_divergence_deg": math.degrees(2 * minimum_theta),
        "effective_beam_quality_m2": theta * math.pi * w0 / wavelength,
        "rayleigh_range_m": z_rayleigh, "nominal_divergence_half_angle_rad": theta,
        "beam_model_status": status,
        "final_beam_diameter_m": None, "final_beam_area_m2": None,
        "beam_area_increase_factor": None, "on_axis_spreading_loss_db": None,
    }
    if valid or allow_invalid:
        diameter = initial_diameter_m * math.hypot(1, distance_m / z_rayleigh)
        area_factor = (diameter / initial_diameter_m)**2
        result.update(final_beam_diameter_m=diameter,
                      final_beam_area_m2=math.pi * (diameter / 2)**2,
                      beam_area_increase_factor=area_factor,
                      on_axis_spreading_loss_db=10 * math.log10(area_factor))
        if not valid:
            result["beam_model_status"] = "unreliable_estimate_" + status
    return result


def divergence_options(config, distance_m):
    """Return regime/angle pairs; overlapping range limits select both presets."""
    mode = config["divergence_mode"]
    if mode == "diffraction_limited":
        return [("diffraction_limited", None)]
    if mode == "specified":
        return [("specified", angle) for angle in config["specified_full_angles_deg"]]
    options = []
    distance_mi = distance_m / MI
    if distance_mi <= config["wide_beam_max_range_mi"]:
        options.extend(("wide", angle) for angle in config["wide_full_angles_deg"])
    if distance_mi >= config["narrow_beam_min_range_mi"]:
        options.extend(("narrow", angle) for angle in config["narrow_full_angles_deg"])
    return options


def endpoints_for(config, start_altitude_ft):
    return ([start_altitude_ft] if config["path_geometry"] == "horizontal"
            else config["slant_end_altitudes_ft"])


def row_count(config, include_impossible=False):
    path_angle_count = sum(
        len(divergence_options(config, range_metres(spec)))
        for start in config["altitudes_ft"] for end in endpoints_for(config, start)
        for spec in config["ranges"]
        if include_impossible or valid_path(range_metres(spec), start, end))
    return (sum(len(v) for v in config["frequencies_ghz"].values())
            * len(config["air_densities_kg_m3"] or config["density_scales"])
            * len(config["weather_cases"]) * path_angle_count
            * len(config["initial_beam_diameters_m"])
            * len(config["pulse_durations_ns"]) * len(config["peak_powers_w"]))


def generate_rows(config, allow_invalid=False):
    """Stream the full Cartesian sweep, reusing gas/cloud calculations."""
    frequencies = [(band, f) for band, values in config["frequencies_ghz"].items() for f in values]
    densities = config["air_densities_kg_m3"] or config["density_scales"]
    for (band, f), altitude, density, weather in itertools.product(
            frequencies, config["altitudes_ft"], densities, config["weather_cases"]):
        for end_altitude in endpoints_for(config, altitude):
            air = integrate_path(altitude, end_altitude, density, weather, config, f,
                                 atmosphere, gas_specific_attenuation, cloud_specific_attenuation)
            dry, wet, cloud = (air["dry_gas_db_per_km"], air["water_vapour_db_per_km"], air["cloud_db_per_km"])
            end_temperature = (config["temperature_k"] if config["temperature_k"] is not None
                               else standard_atmosphere(end_altitude * FT)[0])
            for spec in config["ranges"]:
                distance = range_metres(spec)
                if not valid_path(distance, altitude, end_altitude):
                    continue
                cloud_km = distance / 1000 * air["cloud_path_fraction"]
                gas_loss, cloud_loss = (dry + wet) * distance / 1000, cloud * cloud_km
                loss = gas_loss + cloud_loss
                transmission = 10**(-loss / 10)
                for diameter, (regime, angle) in itertools.product(
                        config["initial_beam_diameters_m"], divergence_options(config, distance)):
                    beam = gaussian_beam(f, diameter, distance, config["beam_quality_m2"],
                                         config["max_paraxial_half_angle_rad"], allow_invalid, angle)
                    for pulse_ns, power in itertools.product(config["pulse_durations_ns"], config["peak_powers_w"]):
                        duration = pulse_ns * 1e-9
                        fractional_offset = 1 / (duration * f * 1e9)
                        flags = ["linear_model_only", "carrier_frequency_only"]
                        flags.append("horizontal_uniform_atmosphere" if config["path_geometry"] == "horizontal"
                                     else "slant_profile_midpoint_integration")
                        if f < 1:
                            flags.append("itu_models_extrapolated_below_1ghz")
                        if fractional_offset > 0.1:
                            flags.append("large_pulse_spectral_offset")
                        if min(air["temperature_k"], end_temperature) < 273.15 and air["liquid_water_g_m3"] > 0:
                            flags.append("supercooled_liquid_assumed_no_ice")
                        if "below_diffraction_limit" in beam["beam_model_status"]:
                            flags.append("requested_divergence_below_diffraction_limit")
                        if beam["nominal_divergence_half_angle_rad"] > config["max_paraxial_half_angle_rad"]:
                            flags.append("beam_outside_paraxial_limit")
                        area = beam["final_beam_area_m2"]
                        yield {
                            "band": band, "frequency_ghz": f, "wavelength_m": C / (f * 1e9),
                            "range_m": distance, "range_ft": distance / FT, "range_mi": distance / MI,
                            "path_geometry": config["path_geometry"],
                            "elevation_deg": elevation_deg(distance, altitude, end_altitude),
                            **air, "cloud_path_km": cloud_km,
                            "gas_attenuation_db": gas_loss, "cloud_attenuation_db": cloud_loss,
                            "atmospheric_attenuation_db": loss, "atmospheric_power_transmission": transmission,
                            **beam, "divergence_regime": regime,
                            "pulse_duration_ns": pulse_ns, "pulse_length_m": C * duration,
                            "first_null_offset_mhz": 1 / duration / 1e6,
                            "first_null_offset_fraction_of_carrier": fractional_offset,
                            "peak_power_w": power, "input_pulse_energy_j": power * duration,
                            "transmitted_peak_power_w": power * transmission,
                            "transmitted_pulse_energy_j": power * duration * transmission,
                            "on_axis_peak_power_density_w_m2": 2 * power * transmission / area if area is not None else None,
                            "model_flags": ";".join(flags),
                        }


def metadata(config, count, allow_invalid):
    data_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (DATA_DIR / "oxygen_lines.csv", DATA_DIR / "water_vapour_lines.csv")}
    return {
        "model_version": MODEL_VERSION, "rows": count, "config": config,
        "rows_skipped_impossible_geometry": row_count(config, True) - count,
        "allow_invalid_beam_estimates": allow_invalid,
        "models": {"gas": "ITU-R P.676-13 Annex 1 eqs. 1–9",
                   "cloud": "ITU-R P.840-9 Annex 1 eqs. 1–10 (local liquid water)",
                   "beam": "Gaussian 1/e^2 intensity diameter, source-plane waist"},
        "coefficient_sha256": data_hashes,
        "sources": [
            "https://www.itu.int/rec/R-REC-P.676-13-202208-I/en",
            "https://www.itu.int/rec/R-REC-P.840-9-202308-I/en",
            "https://github.com/inigodelportillo/ITU-Rpy",
            "https://www.edmundoptics.com/knowledge-center/application-notes/lasers/gaussian-beam-propagation/",
            "https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf",
        ],
        "assumptions": [
            "Straight unobstructed paths; horizontal uses one altitude, slant integrates a linear altitude profile.",
            "Slant density scales follow local ISA; absolute density overrides impose a constant mass density.",
            "Only liquid droplets in cloud/fog; cloud_path_fraction weights occupancy uniformly along the path.",
            "Preset divergence angles are full 1/e^2 opening angles using the paraxial radian slope convention.",
            "Rectangular temporal pulses evaluated at carrier frequency; no pulse distortion or dispersion.",
            "Peak power is launched RF power in the Gaussian mode, not electrical input or EIRP.",
            "Linear propagation; no rain, turbulence, refraction, terrain, curvature, plasma or thermal effects.",
            "No antenna/receiver coupling. Transmitted power is whole-beam power after atmospheric absorption.",
            "Below 1 GHz, ITU attenuation expressions are extrapolations outside their stated scope.",
            "The paraxial angle threshold is a configurable heuristic, not an accuracy guarantee.",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, help="JSON overrides for the built-in sweep")
    parser.add_argument("--output", type=Path, default=Path("atmosphere_results.csv"))
    parser.add_argument("--geometry", choices=["horizontal", "slant"], help="Override path_geometry")
    parser.add_argument("--wide-max-range-mi", type=float, help="Use wide beams up to this distance (inclusive)")
    parser.add_argument("--narrow-min-range-mi", type=float, help="Use narrow beams from this distance (inclusive); may overlap wide beams")
    parser.add_argument("--write-example-config", type=Path, metavar="PATH", help="Write all defaults and exit")
    parser.add_argument("--dry-run", action="store_true", help="Validate and report row count without writing results")
    parser.add_argument("--allow-invalid-beam-estimates", action="store_true", help="Include explicitly flagged Gaussian estimates even when angle constraints fail")
    args = parser.parse_args(argv)
    temporary_paths = []
    try:
        if args.write_example_config:
            args.write_example_config.parent.mkdir(parents=True, exist_ok=True)
            args.write_example_config.write_text(json.dumps(DEFAULT_CONFIG, indent=2) + "\n", encoding="utf-8")
            print(f"Wrote {args.write_example_config}")
            return 0
        config = copy.deepcopy(DEFAULT_CONFIG)
        if args.config:
            overrides = json.loads(args.config.read_text(encoding="utf-8"))
            if not isinstance(overrides, dict):
                raise ValueError("Configuration root must be a JSON object")
            config.update(overrides)
        if args.geometry:
            config["path_geometry"] = args.geometry
        if args.wide_max_range_mi is not None:
            config["wide_beam_max_range_mi"] = args.wide_max_range_mi
        if args.narrow_min_range_mi is not None:
            config["narrow_beam_min_range_mi"] = args.narrow_min_range_mi
        validate_config(config)
        count = row_count(config)
        if args.dry_run:
            print(f"Validated configuration: {count:,} rows (no files written).")
            return 0
        meta = metadata(config, count, args.allow_invalid_beam_estimates)
        # Detect missing/malformed tables before preparing output files.
        spectral_lines("oxygen_lines.csv")
        spectral_lines("water_vapour_lines.csv")
        sidecar = args.output.with_suffix(args.output.suffix + ".metadata.json")
        protected = {Path(__file__).resolve(), Path(__file__).resolve().parent / "rf_geometry.py", *(DATA_DIR / name for name in
                      ["oxygen_lines.csv", "water_vapour_lines.csv", "ITU-Rpy-LICENSE.txt"])}
        if args.config:
            protected.add(args.config.resolve())
        if args.output.resolve() in {p.resolve() for p in protected} or sidecar.resolve() in {p.resolve() for p in protected}:
            raise ValueError("Output would overwrite an input, script or coefficient file")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        invalid = 0
        below_diffraction = 0
        outside_paraxial = 0
        status_counts = {}
        with tempfile.NamedTemporaryFile(mode="w", newline="", encoding="utf-8",
                                         prefix=".rf-atmosphere-", suffix=".csv.tmp",
                                         dir=args.output.parent, delete=False) as handle:
            csv_temporary = Path(handle.name)
            temporary_paths.append(csv_temporary)
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for row in generate_rows(config, args.allow_invalid_beam_estimates):
                if any(isinstance(v, (float, int)) and not math.isfinite(v) for v in row.values()):
                    raise ValueError("Inputs produced nonfinite results; reduce extreme density, water content, beam quality or range")
                invalid += row["beam_model_status"] != "paraxial_estimate"
                below_diffraction += "requested_divergence_below_diffraction_limit" in row["model_flags"]
                outside_paraxial += "beam_outside_paraxial_limit" in row["model_flags"]
                status = row["beam_model_status"]
                status_counts[status] = status_counts.get(status, 0) + 1
                # Keep small attenuation values visible, with useful numeric precision.
                writer.writerow({k: format(v, ".12g") if isinstance(v, float) else v for k, v in row.items()})
        meta["rows_outside_paraxial_limit"] = outside_paraxial
        meta["rows_below_diffraction_limit"] = below_diffraction
        meta["rows_without_valid_beam_prediction"] = invalid
        meta["beam_model_status_counts"] = status_counts
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", prefix=".rf-atmosphere-",
                                         suffix=".json.tmp", dir=sidecar.parent, delete=False) as handle:
            metadata_temporary = Path(handle.name)
            temporary_paths.append(metadata_temporary)
            handle.write(json.dumps(meta, indent=2, allow_nan=False) + "\n")
        csv_temporary.replace(args.output)
        metadata_temporary.replace(sidecar)
        print(f"Wrote {count:,} rows to {args.output}")
        print(f"Model settings and sources: {sidecar}")
        if invalid:
            action = "flagged estimates included" if args.allow_invalid_beam_estimates else "beam predictions left blank"
            print(f"{invalid:,} rows fail Gaussian angle constraints: {action}.")
        if meta["rows_skipped_impossible_geometry"]:
            print(f"Skipped {meta['rows_skipped_impossible_geometry']:,} rows with impossible slant geometry.")
        return 0
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as error:
        parser.exit(2, f"Error: {error}\n")
    finally:
        for path in temporary_paths:
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
