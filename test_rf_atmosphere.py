"""Physical invariants, published gas reference, and end-to-end CSV checks."""

import copy
import csv
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import rf_atmosphere as rf


def small_config():
    config = copy.deepcopy(rf.DEFAULT_CONFIG)
    config.update(frequencies_ghz={"X": [10]}, altitudes_ft=[0], density_scales=[1],
                  divergence_mode="diffraction_limited",
                  weather_cases=[config["weather_cases"][0]],
                  ranges=[{"value": 1, "unit": "mi"}],
                  initial_beam_diameters_m=[1], pulse_durations_ns=[10, 1000],
                  peak_powers_w=[1e3, 1e9])
    return config


class PropagationTests(unittest.TestCase):
    def test_published_itu_gas_reference(self):
        # ITU-Rpy validation of P.676-12 (same exact formulas/tables as -13).
        # https://itu-rpy.readthedocs.io/en/latest/validation/itur676_12.html
        dry, wet = rf.gas_specific_attenuation(22, 1013.25, 7.5, 288.15)
        self.assertAlmostEqual(dry + wet, 0.187337, delta=0.0000005)

    def test_no_water_vapour_has_no_wet_absorption(self):
        dry, wet = rf.gas_specific_attenuation(22.23508, 1013.25, 0, 288.15)
        self.assertGreater(dry, 0)
        self.assertEqual(wet, 0)

    def test_cloud_zero_and_linear_liquid_water_scaling(self):
        self.assertEqual(rf.cloud_specific_attenuation(35, 288.15, 0), 0)
        one = rf.cloud_specific_attenuation(35, 288.15, 0.25)
        two = rf.cloud_specific_attenuation(35, 288.15, 0.5)
        self.assertGreater(one, 0)
        self.assertAlmostEqual(two, 2 * one)

    def test_rayleigh_distance_and_far_field(self):
        beam = rf.gaussian_beam(10, 1, 0)
        rayleigh = beam["rayleigh_range_m"]
        at_rayleigh = rf.gaussian_beam(10, 1, rayleigh)
        self.assertAlmostEqual(at_rayleigh["final_beam_diameter_m"], math.sqrt(2))
        far = rf.gaussian_beam(10, 1, rayleigh * 10000)
        expected = 2 * beam["nominal_divergence_half_angle_rad"] * rayleigh * 10000
        self.assertAlmostEqual(far["final_beam_diameter_m"] / expected, 1, places=7)

    def test_invalid_uhf_beam_is_blank_unless_explicitly_requested(self):
        beam = rf.gaussian_beam(0.5, 0.2, 1000)
        self.assertEqual(beam["beam_model_status"], "outside_paraxial_limit")
        self.assertIsNone(beam["final_beam_diameter_m"])
        override = rf.gaussian_beam(0.5, 0.2, 1000, allow_invalid=True)
        self.assertGreater(override["final_beam_diameter_m"], 0)
        self.assertIn("unreliable", override["beam_model_status"])

    def test_density_is_mass_density_and_pressure_is_consistent(self):
        config = small_config()
        weather = config["weather_cases"][0]
        base = rf.atmosphere(0, 1, weather, config)
        scaled = rf.atmosphere(0, 0.8, weather, config)
        self.assertAlmostEqual(scaled["air_density_kg_m3"], 0.8 * base["air_density_kg_m3"])
        e = base["water_vapour_density_g_m3"] * base["temperature_k"] / 216.7
        self.assertAlmostEqual(base["dry_pressure_hpa"] + e, base["total_pressure_hpa"])
        self.assertAlmostEqual(base["total_pressure_hpa"], 1013.25)
        config.update(density_scales=None, air_densities_kg_m3=[0.9])
        explicit = rf.atmosphere(0, 0.9, weather, config)
        self.assertAlmostEqual(explicit["air_density_kg_m3"], 0.9)

    def test_altitude_reduces_reference_density(self):
        config = small_config()
        weather = config["weather_cases"][0]
        low = rf.atmosphere(0, 1, weather, config)
        high = rf.atmosphere(10000, 1, weather, config)
        self.assertLess(high["air_density_kg_m3"], low["air_density_kg_m3"])
        self.assertLess(high["total_pressure_hpa"], low["total_pressure_hpa"])

    def test_zero_path_and_path_length_scaling(self):
        config = small_config()
        config["weather_cases"] = [rf.DEFAULT_CONFIG["weather_cases"][2]]
        config["ranges"] = [{"value": v, "unit": "km"} for v in [0, 1, 2]]
        config["pulse_durations_ns"] = [10]
        config["peak_powers_w"] = [1e3]
        zero, one, two = list(rf.generate_rows(config))
        self.assertEqual(zero["atmospheric_attenuation_db"], 0)
        self.assertEqual(zero["atmospheric_power_transmission"], 1)
        self.assertEqual(zero["final_beam_diameter_m"], 1)
        self.assertAlmostEqual(two["gas_attenuation_db"], 2 * one["gas_attenuation_db"])
        self.assertAlmostEqual(two["cloud_attenuation_db"], 2 * one["cloud_attenuation_db"])
        config["weather_cases"][0] = dict(config["weather_cases"][0], cloud_path_fraction=0.5)
        half = list(rf.generate_rows(config))[1]
        self.assertAlmostEqual(half["cloud_attenuation_db"], 0.5 * one["cloud_attenuation_db"])
        self.assertEqual(half["gas_attenuation_db"], one["gas_attenuation_db"])

    def test_energy_power_and_gaussian_intensity(self):
        config = small_config()
        rows = list(rf.generate_rows(config))
        for row in rows:
            transmission = row["atmospheric_power_transmission"]
            self.assertAlmostEqual(row["transmitted_peak_power_w"] / row["peak_power_w"], transmission)
            self.assertAlmostEqual(row["transmitted_pulse_energy_j"] / row["input_pulse_energy_j"], transmission)
            self.assertAlmostEqual(row["on_axis_peak_power_density_w_m2"] * row["final_beam_area_m2"] / 2,
                                   row["transmitted_peak_power_w"])
        self.assertEqual(len({r["final_beam_diameter_m"] for r in rows}), 1)
        self.assertEqual(len({r["atmospheric_attenuation_db"] for r in rows}), 1)
        self.assertAlmostEqual(rows[0]["pulse_length_m"], 2.99792458)
        self.assertAlmostEqual(rows[0]["first_null_offset_mhz"], 100)
        self.assertAlmostEqual(rows[-1]["input_pulse_energy_j"], 1000)

    def test_default_sweep_count_and_parameter_bounds(self):
        config = copy.deepcopy(rf.DEFAULT_CONFIG)
        rf.validate_config(config)
        self.assertEqual(rf.row_count(config), 122472)
        for key, invalid in [("peak_powers_w", [float("nan")]),
                             ("pulse_durations_ns", [9]),
                             ("initial_beam_diameters_m", [0]),
                             ("altitudes_ft", [40000]),
                             ("density_scales", [0])]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                modified = dict(config, **{key: invalid})
                rf.validate_config(modified)

    def test_csv_export_and_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            config_path = Path(folder) / "config.json"
            output = Path(folder) / "results.csv"
            config_path.write_text(json.dumps(small_config()))
            self.assertEqual(rf.main(["--config", str(config_path), "--output", str(output)]), 0)
            with output.open(newline="") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(reader.fieldnames, rf.FIELDS)
                rows = list(reader)
            self.assertEqual(len(rows), 4)
            meta = json.loads(output.with_suffix(".csv.metadata.json").read_text())
            self.assertEqual(meta["rows"], 4)
            self.assertEqual(meta["rows_outside_paraxial_limit"], 0)
            self.assertEqual(len(meta["coefficient_sha256"]), 2)
            self.assertLess(float(rows[0]["transmitted_peak_power_w"]), 1000)
            self.assertAlmostEqual(float(rows[0]["range_ft"]), 5280)

    def test_uhf_model_and_pulse_flags(self):
        config = small_config()
        config["frequencies_ghz"] = {"UHF": [0.5]}
        config["initial_beam_diameters_m"] = [0.2]
        row = next(rf.generate_rows(config))
        self.assertIn("extrapolated_below_1ghz", row["model_flags"])
        self.assertIn("large_pulse_spectral_offset", row["model_flags"])
        self.assertIsNone(row["on_axis_peak_power_density_w_m2"])

    def test_range_overflow_rejected(self):
        with self.assertRaises(ValueError):
            rf.range_metres({"value": 1e308, "unit": "km"})

    def test_generation_failure_preserves_previous_output(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "results.csv"
            output.write_text("previous results\n")
            sidecar = output.with_suffix(".csv.metadata.json")
            sidecar.write_text("previous metadata\n")
            with patch.object(rf, "generate_rows", side_effect=ValueError("test calculation failure")):
                with self.assertRaises(SystemExit):
                    rf.main(["--output", str(output)])
            self.assertEqual(output.read_text(), "previous results\n")
            self.assertEqual(sidecar.read_text(), "previous metadata\n")
            self.assertEqual(sorted(p.name for p in Path(folder).iterdir()),
                             ["results.csv", "results.csv.metadata.json"])

    def test_overlapping_beam_ranges_include_both_at_inclusive_limits(self):
        config = copy.deepcopy(rf.DEFAULT_CONFIG)
        for miles, expected in [(0.5, 3), (1, 6), (3, 6), (5, 6), (6, 3)]:
            with self.subTest(miles=miles):
                options = rf.divergence_options(config, miles * rf.MI)
                self.assertEqual(len(options), expected)
        config.update(wide_beam_max_range_mi=10, narrow_beam_min_range_mi=2)
        self.assertEqual(len(rf.divergence_options(config, 3 * rf.MI)), 6)
        self.assertEqual(len(rf.divergence_options(config, rf.MI)), 3)

    def test_requested_divergence_respects_diffraction_floor(self):
        invalid = rf.gaussian_beam(0.5, 2, 10000, full_angle_deg=5)
        self.assertEqual(invalid["beam_model_status"], "below_diffraction_limit")
        self.assertIsNone(invalid["final_beam_diameter_m"])
        self.assertGreater(invalid["diffraction_min_full_divergence_deg"], 21)
        valid = rf.gaussian_beam(10, 1, 10000, full_angle_deg=3)
        self.assertEqual(valid["beam_model_status"], "paraxial_estimate")
        self.assertGreater(valid["effective_beam_quality_m2"], 1)
        expected_diameter = 2 * math.hypot(0.5, 10000 * math.radians(1.5))
        self.assertAlmostEqual(valid["final_beam_diameter_m"], expected_diameter)
        other_frequency = rf.gaussian_beam(35, 1, 10000, full_angle_deg=3)
        self.assertEqual(valid["final_beam_diameter_m"], other_frequency["final_beam_diameter_m"])

    def test_slant_counts_and_impossible_ranges(self):
        config = copy.deepcopy(rf.DEFAULT_CONFIG)
        config["path_geometry"] = "slant"
        rf.validate_config(config)
        self.assertEqual(rf.row_count(config), 193914)
        self.assertEqual(rf.row_count(config, True) - rf.row_count(config), 51030)
        tiny = small_config()
        tiny.update(path_geometry="slant", slant_end_altitudes_ft=[10000],
                    ranges=[{"value": 1000, "unit": "ft"}, {"value": 50, "unit": "mi"}])
        rows = list(rf.generate_rows(tiny))
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(r["end_altitude_ft"] == 10000 for r in rows))
        self.assertTrue(all(r["range_mi"] == 50 for r in rows))

    def test_slant_integrated_loss_reversal_and_convergence(self):
        config = small_config()
        config.update(path_geometry="slant", slant_end_altitudes_ft=[10000],
                      ranges=[{"value": 50, "unit": "mi"}])
        forward = next(rf.generate_rows(config))
        self.assertLess(forward["path_mean_temperature_k"], forward["temperature_k"])
        self.assertGreater(forward["elevation_deg"], 0)
        config.update(altitudes_ft=[10000], slant_end_altitudes_ft=[0])
        reverse = next(rf.generate_rows(config))
        self.assertAlmostEqual(forward["atmospheric_attenuation_db"], reverse["atmospheric_attenuation_db"])
        config["slant_integration_steps"] = 256
        fine = next(rf.generate_rows(config))
        self.assertAlmostEqual(reverse["atmospheric_attenuation_db"] / fine["atmospheric_attenuation_db"], 1, places=4)


if __name__ == "__main__":
    unittest.main()
