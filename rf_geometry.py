"""Straight horizontal/slant path helpers, independent of atmospheric models.

The slant path follows a linear altitude profile. Its midpoint quadrature
returns length-averaged atmospheric properties and attenuation rates. Multiply
the rates by the actual path length elsewhere to obtain total attenuation.
"""

import math

FT = 0.3048


def valid_path(range_m, start_altitude_ft, end_altitude_ft):
    """Whether a straight path can span the endpoint height difference.

    Range is the total path length, not ground distance. A tiny tolerance
    permits floating-point roundoff in a nominally vertical path. A zero-length
    path requires exactly equal endpoint altitudes.
    """
    values = (range_m, start_altitude_ft, end_altitude_ft)
    if any(not math.isfinite(value) for value in values) or range_m < 0:
        return False
    if range_m == 0:
        return start_altitude_ft == end_altitude_ft
    height_difference_m = abs(end_altitude_ft - start_altitude_ft) * FT
    return (range_m >= height_difference_m
            or math.isclose(range_m, height_difference_m,
                            rel_tol=1e-12, abs_tol=1e-9))


def elevation_deg(range_m, start_altitude_ft, end_altitude_ft):
    """Signed straight-path elevation in degrees, or 0 for a zero-length path."""
    if not valid_path(range_m, start_altitude_ft, end_altitude_ft):
        raise ValueError("Path range must span the endpoint altitude difference")
    if range_m == 0:
        return 0.0
    ratio = (end_altitude_ft - start_altitude_ft) * FT / range_m
    return math.degrees(math.asin(max(-1.0, min(1.0, ratio))))


def integrate_path(start_altitude_ft, end_altitude_ft, density, weather, config,
                   f_ghz, atmosphere_fn, gas_fn, cloud_fn):
    """Return source properties and path-mean loss rates for uniform sampling.

    Callbacks use these signatures:
      atmosphere_fn(altitude_ft, density, weather, config) -> atmosphere dict
      gas_fn(f_ghz, dry_pressure_hpa, vapour_g_m3, temperature_k) -> (dry, wet)
      cloud_fn(f_ghz, temperature_k, liquid_water_g_m3) -> cloud loss rate

    Loss rates are in dB/km. The liquid-water rate is the in-cloud rate;
    cloud-path occupancy is applied by the caller. For unequal endpoint
    altitudes, midpoint quadrature uses ``slant_integration_steps`` (default64).
    """
    source = atmosphere_fn(start_altitude_ft, density, weather, config)

    def local_rates(air):
        dry, wet = gas_fn(f_ghz, air["dry_pressure_hpa"],
                          air["water_vapour_density_g_m3"], air["temperature_k"])
        cloud = cloud_fn(f_ghz, air["temperature_k"], air["liquid_water_g_m3"])
        return dry, wet, cloud

    if start_altitude_ft == end_altitude_ft:
        mean_temperature = source["temperature_k"]
        mean_density = source["air_density_kg_m3"]
        dry, wet, cloud = local_rates(source)
    else:
        steps = config.get("slant_integration_steps", 64)
        if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
            raise ValueError("slant_integration_steps must be a positive integer")
        samples = []
        altitude_span = end_altitude_ft - start_altitude_ft
        for i in range(steps):
            altitude = start_altitude_ft + (i + 0.5) * altitude_span / steps
            air = atmosphere_fn(altitude, density, weather, config)
            samples.append((*local_rates(air), air["temperature_k"],
                            air["air_density_kg_m3"]))
        dry, wet, cloud, mean_temperature, mean_density = (
            math.fsum(values) / steps for values in zip(*samples))

    return {
        **source,
        "end_altitude_ft": end_altitude_ft,
        "end_altitude_m": end_altitude_ft * FT,
        "path_mean_temperature_k": mean_temperature,
        "path_mean_air_density_kg_m3": mean_density,
        "dry_gas_db_per_km": dry,
        "water_vapour_db_per_km": wet,
        "cloud_db_per_km": cloud,
    }
