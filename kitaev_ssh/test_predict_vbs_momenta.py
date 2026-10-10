#!/usr/bin/env python3
"""Focused generator/unit tests; the C++ parser bridge lives in the e2e test."""
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from predict_vbs_momenta import (cartesian_distance, momentum_distance,
                                 momentum_separation, nearest_grid,
                                 nearest_grid_set,
                                 reciprocal_cartesian, uniform_background_solution)


HERE = Path(__file__).resolve().parent
PREDICTOR = HERE / "predict_vbs_momenta.py"
GEOMETRY = HERE / "compare_sheared_geometries.py"


class PredictorTests(unittest.TestCase):
    def run_predictor(self, *arguments):
        directory = tempfile.TemporaryDirectory(prefix="predict-vbs-test-")
        command = [sys.executable, str(PREDICTOR), *map(str, arguments),
                   "--output-dir", directory.name, "--prefix", "case"]
        environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
        completed = subprocess.run(command, text=True, capture_output=True, env=environment)
        self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
        report = json.loads((Path(directory.name) / "case_report.json").read_text())
        qfile = (Path(directory.name) / "case_qfile.txt").read_bytes()
        return directory, report, qfile

    def test_uniform_background_converges_with_project_conventions(self):
        result = uniform_background_solution((1., 1., .5), .3, 1.7 * .8**2,
                                             tolerance=1e-9, integration_grid=96)
        self.assertTrue(result["converged"])
        self.assertLessEqual(result["final_residual_max_abs_X"], 1e-9)
        self.assertAlmostEqual(result["Xbar"][0], result["Xbar"][1], places=13)
        self.assertTrue(all(x > 0 for x in result["Xbar"]))
        self.assertTrue(all(e > j for e, j in zip(result["Jeff"], (1., 1., .5))))

    def test_cartesian_distance_uses_actual_reciprocal_basis_and_collinear_additivity(self):
        bare = (0.5804306232551663, 0.4195693767448336)
        uniform = (0.5806372615976502, 0.4193627384023497)
        grid = (7. / 12., 5. / 12.)
        bare_uniform = momentum_distance(bare, uniform)
        uniform_grid = momentum_distance(uniform, grid)
        bare_grid = momentum_distance(bare, grid)
        self.assertAlmostEqual(bare_uniform["cartesian_distance"] +
                               uniform_grid["cartesian_distance"],
                               bare_grid["cartesian_distance"], places=14)
        self.assertAlmostEqual(bare_uniform["cartesian_distance"],
                               0.001499201976894925, places=14)
        direct = reciprocal_cartesian(bare[0] - uniform[0], bare[1] - uniform[1])
        self.assertAlmostEqual(math.hypot(*direct), bare_uniform["cartesian_distance"], places=15)
        self.assertAlmostEqual(cartesian_distance(*bare_uniform["fractional_delta"]),
                               bare_uniform["cartesian_distance"], places=15)

    def test_periodic_and_minus_momentum_distance(self):
        q = (.2, .3)
        self.assertLess(momentum_distance(q, (1.2, -.7))["cartesian_distance"], 1e-14)
        self.assertLess(momentum_distance((.5, .5), (-.5, -.5))["cartesian_distance"], 1e-14)
        self.assertGreater(momentum_distance(q, (-q[0], -q[1]))["cartesian_distance"], 1e-6)
        self.assertAlmostEqual(momentum_distance(q, (-q[0], -q[1]))["cartesian_distance"],
                               momentum_distance((-q[0], -q[1]), q)["cartesian_distance"], places=14)

    def test_nearest_ties_are_complete_stable_and_not_qfile_deduplication(self):
        bare = (.5804306232551663, .41956937674483363)
        uniform = (.5806552317560235, .4193447682439766)
        bare_tie = nearest_grid_set(bare, 6, 6, 0)
        uniform_tie = nearest_grid_set(uniform, 6, 6, 0)
        self.assertEqual([p["index"] for p in bare_tie], [20, 27])
        self.assertEqual([p["index"] for p in uniform_tie], [20, 27])
        self.assertAlmostEqual(bare_tie[0]["cartesian_distance"], bare_tie[1]["cartesian_distance"], places=14)
        broken = nearest_grid_set((uniform[0] + 1e-4, uniform[1]), 6, 6, 0)
        self.assertEqual(len(broken), 1)
        periodic = nearest_grid_set((uniform[0] + 1., uniform[1] - 1.), 6, 6, 0)
        self.assertEqual([p["index"] for p in periodic], [20, 27])

    def test_geometry_report_reuses_target_distance_definition(self):
        with tempfile.TemporaryDirectory(prefix="geometry-distance-test-") as directory:
            command = [sys.executable, str(GEOMETRY), "--uniform-effective", "--Jx", "1",
                       "--Jy", "1", "--Jz", ".5", "--g", ".3", "--mass", "1.7",
                       "--omega-ph", ".8", "--uniform-grid", "400", "--geometry", "12,3,3",
                       "--output-dir", directory, "--prefix", "case"]
            completed = subprocess.run(command, text=True, capture_output=True,
                                       env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            report = json.loads((Path(directory) / "case_report.json").read_text())
            comparison = report["target_comparisons"][0]
            self.assertAlmostEqual(comparison["cartesian_distance"],
                                   0.001499201976894925, places=14)
            row = report["reports"][0]["targets"][0]
            self.assertAlmostEqual(row["cartesian_mismatch"],
                                   0.021059734777342673, places=14)

    def test_geometry_report_keeps_symmetric_nearest_ties(self):
        with tempfile.TemporaryDirectory(prefix="geometry-tie-test-") as directory:
            command = [sys.executable, str(GEOMETRY), "--uniform-effective", "--Jx", "1",
                       "--Jy", "1", "--Jz", ".5", "--g", ".3", "--mass", "1",
                       "--omega-ph", "1", "--uniform-grid", "400", "--geometry", "6,6,0",
                       "--output-dir", directory, "--prefix", "case"]
            completed = subprocess.run(command, text=True, capture_output=True,
                                       env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            report = json.loads((Path(directory) / "case_report.json").read_text())
            targets = {row["label"]: row for row in report["reports"][0]["targets"]}
            for label in ("bare_Qp", "uniform_effective_Qp"):
                self.assertEqual([p["index"] for p in targets[label]["nearest_candidates"]], [20, 27])

    def test_actual_coordinates_deduplicate_but_nearest_labels_do_not(self):
        directory, report, _ = self.run_predictor(
            "--Lx", 2, "--Ly", 2, "--Jz", .05, "--neighbors", 1,
            "--extra-q", "left,0.01,0.01", "--extra-q", "right,-0.01,-0.01")
        self.addCleanup(directory.cleanup)
        points = report["q_file_measurement_points"]
        shared = [p for p in points if "bare_Qp_g3" in p["aliases"]]
        self.assertEqual(len(shared), 1)
        self.assertIn("bare_Qm_g3", shared[0]["aliases"])
        left = next(p for p in points if "user_left" in p["aliases"])
        right = next(p for p in points if "user_right" in p["aliases"])
        self.assertNotEqual(left["index"], right["index"])
        self.assertGreater(momentum_separation(left["reciprocal_fraction"],
                                               right["reciprocal_fraction"]), 1e-3)
        self.assertEqual(nearest_grid(left["reciprocal_fraction"], 2, 2, 0, 1)[0]["index"], 0)
        self.assertEqual(nearest_grid(right["reciprocal_fraction"], 2, 2, 0, 1)[0]["index"], 0)

    def test_kekule_uniform_and_special_minus_aliases(self):
        directory, report, _ = self.run_predictor(
            "--Lx", 6, "--Ly", 6, "--Jx", 1, "--Jy", 1, "--Jz", 1,
            "--neighbors", 1, "--uniform-effective", "--g", 0,
            "--spring-k", 1, "--uniform-grid", 64,
            "--extra-q", "special_Q,.5,.5",
            "--extra-q", "special_minus,-.5,-.5")
        self.addCleanup(directory.cleanup)
        points = report["q_file_measurement_points"]
        self.assertTrue(any(any(a.startswith("bare_Q") for a in p["aliases"])
                            and any(a.startswith("kek_K") for a in p["aliases"])
                            for p in points))
        self.assertTrue(any(any(a.startswith("bare_Q") for a in p["aliases"])
                            and any(a.startswith("uniform_effective_Q") for a in p["aliases"])
                            for p in points))
        special = next(p for p in points if "user_special_Q" in p["aliases"])
        self.assertIn("user_special_minus", special["aliases"])
        self.assertEqual(special["reciprocal_fraction"], [.5, .5])

    def test_extra_q_grid_metadata_comes_from_geometry(self):
        directory, report, _ = self.run_predictor(
            "--Lx", 12, "--Ly", 3, "--shear", 3,
            "--extra-q", "grid,0.4166666666666667,0.5833333333333333",
            "--extra-q", "off,.137,.271", "--neighbors", 1)
        self.addCleanup(directory.cleanup)
        candidates = {row["label"]: row for row in report["candidates"]}
        self.assertTrue(candidates["user_grid"]["on_grid"])
        self.assertEqual((candidates["user_grid"]["n1"], candidates["user_grid"]["n2"]), (5, 0))
        self.assertFalse(candidates["user_off"]["on_grid"])

    def test_qfile_order_is_deterministic(self):
        args = ("--Lx", 12, "--Ly", 3, "--shear", 3, "--neighbors", 2,
                "--uniform-effective", "--g", .3, "--mass", 1.7,
                "--omega-ph", .8, "--uniform-grid", 64, "--extra-q", "zero,0,0")
        first, first_report, first_qfile = self.run_predictor(*args)
        second, second_report, second_qfile = self.run_predictor(*args)
        self.addCleanup(first.cleanup)
        self.addCleanup(second.cleanup)
        self.assertEqual(first_qfile, second_qfile)
        self.assertEqual(first_report["q_file_measurement_points"],
                         second_report["q_file_measurement_points"])

    def test_nonconvergence_is_an_error(self):
        with tempfile.TemporaryDirectory(prefix="predict-vbs-fail-") as directory:
            command = [sys.executable, str(PREDICTOR), "--Lx", "6", "--Ly", "6",
                       "--uniform-effective", "--g", ".3", "--spring-k", "1.088",
                       "--uniform-max-iterations", "1", "--uniform-tolerance", "1e-16",
                       "--uniform-grid", "64", "--output-dir", directory, "--prefix", "fail"]
            completed = subprocess.run(command, text=True, capture_output=True,
                                       env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("did not converge", completed.stderr + completed.stdout)
            self.assertFalse((Path(directory) / "fail_qfile.txt").exists())


if __name__ == "__main__":
    unittest.main()
