"""Geometry and path quadrature checks with analytically integrable models."""

import math
import unittest

import rf_geometry as geometry


def linear_atmosphere(altitude_ft, density, weather, config):
    height = altitude_ft * geometry.FT
    return {
        "altitude_ft": altitude_ft,
        "altitude_m": height,
        "temperature_k": 280 + 0.001 * height,
        "air_density_kg_m3": density + 0.0001 * height,
        "dry_pressure_hpa": 900 + 0.002 * height,
        "water_vapour_density_g_m3": 3 + 0.003 * height,
        "liquid_water_g_m3": weather["liquid_water_g_m3"],
    }


def linear_gas(f_ghz, dry_pressure_hpa, vapour_g_m3, temperature_k):
    return dry_pressure_hpa / 100, vapour_g_m3 / 10


def linear_cloud(f_ghz, temperature_k, liquid_water_g_m3):
    return temperature_k * liquid_water_g_m3 / 1000


class GeometryTests(unittest.TestCase):
    def test_range_is_slant_length_and_must_span_height(self):
        height = 1000 * geometry.FT
        self.assertFalse(geometry.valid_path(height - 1, 0, 1000))
        self.assertTrue(geometry.valid_path(height, 0, 1000))
        self.assertTrue(geometry.valid_path(height + 1, 0, 1000))
        self.assertTrue(geometry.valid_path(height, 1000, 0))

    def test_invalid_ranges_and_zero_length(self):
        for distance in (-1, math.inf, math.nan):
            self.assertFalse(geometry.valid_path(distance, 0, 0))
        self.assertFalse(geometry.valid_path(1000, 0, math.inf))
        self.assertTrue(geometry.valid_path(0, 1000, 1000))
        self.assertFalse(geometry.valid_path(0, 1000, 1000.000000001))

    def test_elevation_and_roundoff_at_vertical_path(self):
        self.assertAlmostEqual(geometry.elevation_deg(1000, 0, 500 / geometry.FT), 30)
        self.assertAlmostEqual(geometry.elevation_deg(1000, 500 / geometry.FT, 0), -30)
        self.assertEqual(geometry.elevation_deg(0, 1000, 1000), 0)
        self.assertEqual(geometry.elevation_deg(geometry.FT * 1000, 0, 1000), 90)
        self.assertEqual(geometry.elevation_deg(geometry.FT * 1000 - 1e-10, 0, 1000), 90)
        with self.assertRaises(ValueError):
            geometry.elevation_deg(100, 0, 1000)


class IntegrationTests(unittest.TestCase):
    weather = {"liquid_water_g_m3": 0.5}

    def integrate(self, start, end, config=None, atmosphere_fn=linear_atmosphere):
        return geometry.integrate_path(start, end, 1.2, self.weather, config or {},
                                       10, atmosphere_fn, linear_gas, linear_cloud)

    def test_horizontal_matches_local_atmosphere(self):
        source = linear_atmosphere(1000, 1.2, self.weather, {})
        result = self.integrate(1000, 1000)
        for key, value in source.items():
            self.assertEqual(result[key], value)
        self.assertEqual(result["end_altitude_ft"], 1000)
        self.assertEqual(result["path_mean_temperature_k"], source["temperature_k"])
        self.assertEqual(result["path_mean_air_density_kg_m3"], source["air_density_kg_m3"])
        self.assertEqual(result["dry_gas_db_per_km"], source["dry_pressure_hpa"] / 100)

    def test_midpoints_integrate_linear_profiles_exactly(self):
        midpoint = linear_atmosphere(6000, 1.2, self.weather, {})
        for steps in (1, 2, 64):
            result = self.integrate(2000, 10000, {"slant_integration_steps": steps})
            self.assertAlmostEqual(result["path_mean_temperature_k"], midpoint["temperature_k"])
            self.assertAlmostEqual(result["path_mean_air_density_kg_m3"], midpoint["air_density_kg_m3"])
            self.assertAlmostEqual(result["dry_gas_db_per_km"], midpoint["dry_pressure_hpa"] / 100)
            self.assertAlmostEqual(result["water_vapour_db_per_km"], midpoint["water_vapour_density_g_m3"] / 10)
            self.assertAlmostEqual(result["cloud_db_per_km"], midpoint["temperature_k"] * 0.5 / 1000)
            self.assertEqual(result["altitude_ft"], 2000)
            self.assertEqual(result["end_altitude_ft"], 10000)

    def test_reversing_path_preserves_means_and_changes_source(self):
        upward = self.integrate(2000, 10000)
        downward = self.integrate(10000, 2000)
        keys = ("path_mean_temperature_k", "path_mean_air_density_kg_m3",
                "dry_gas_db_per_km", "water_vapour_db_per_km", "cloud_db_per_km")
        for key in keys:
            self.assertAlmostEqual(upward[key], downward[key])
        self.assertNotEqual(upward["altitude_ft"], downward["altitude_ft"])

    def test_midpoint_convergence_for_quadratic_profile(self):
        def quadratic_atmosphere(altitude_ft, density, weather, config):
            result = linear_atmosphere(altitude_ft, density, weather, config)
            result["air_density_kg_m3"] = (altitude_ft / 1000)**2
            return result

        coarse = self.integrate(0, 1000, {"slant_integration_steps": 2}, quadratic_atmosphere)
        fine = self.integrate(0, 1000, {"slant_integration_steps": 4}, quadratic_atmosphere)
        exact = 1 / 3
        coarse_error = exact - coarse["path_mean_air_density_kg_m3"]
        fine_error = exact - fine["path_mean_air_density_kg_m3"]
        self.assertAlmostEqual(coarse_error / fine_error, 4)

    def test_invalid_step_counts(self):
        for steps in (0, -1, True, 1.5, "64"):
            with self.assertRaises(ValueError):
                self.integrate(0, 1000, {"slant_integration_steps": steps})


if __name__ == "__main__":
    unittest.main()
