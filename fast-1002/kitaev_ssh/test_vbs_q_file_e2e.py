#!/usr/bin/env python3
"""predict_vbs_momenta.py -> q-file -> production C++ parser metadata checks."""
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
PREDICTOR = HERE / "predict_vbs_momenta.py"
VALIDATOR = os.environ.get("VBS_Q_VALIDATE_BIN")


class QFileEndToEndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not VALIDATOR:
            raise unittest.SkipTest("set VBS_Q_VALIDATE_BIN to the built kitaev_ssh_vbs_q_validate")
        cls.validator = Path(VALIDATOR).resolve()
        if not cls.validator.is_file():
            raise RuntimeError(f"validator binary does not exist: {cls.validator}")

    def generate_and_parse(self, name, Lx, Ly, shear, *extra):
        directory = tempfile.TemporaryDirectory(prefix=f"vbs-q-e2e-{name}-")
        self.addCleanup(directory.cleanup)
        prefix = Path(directory.name) / name
        command = [sys.executable, str(PREDICTOR), "--Lx", str(Lx), "--Ly", str(Ly),
                   "--shear", str(shear), "--output-dir", directory.name, "--prefix", name,
                   *map(str, extra)]
        generated = subprocess.run(command, text=True, capture_output=True,
                                   env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        self.assertEqual(generated.returncode, 0, generated.stderr + generated.stdout)
        report = json.loads(prefix.with_name(name + "_report.json").read_text())
        qfile = prefix.with_name(name + "_qfile.txt")
        parsed = subprocess.run([str(self.validator), "--inspect-q-file", str(Lx), str(Ly),
                                 str(shear), str(qfile)], text=True, capture_output=True)
        self.assertEqual(parsed.returncode, 0, parsed.stderr + parsed.stdout)
        cpp = json.loads(parsed.stdout)["points"]
        expected = report["q_file_measurement_points"]
        self.assertEqual(len(cpp), len(expected))
        for actual, wanted in zip(cpp, expected):
            self.assertEqual(actual["index"], wanted["index"])
            self.assertEqual(actual["label"], wanted["label"])
            self.assertEqual(actual["on_grid"], wanted["on_grid"])
            self.assertEqual(actual["n"], wanted["n"])
            for key in ("reciprocal_fraction", "phase_fraction"):
                self.assertLess(math.dist(actual[key], wanted[key]), 2e-14)
        return report, cpp

    def test_required_geometries_and_metadata(self):
        common = ("--Jx", 1, "--Jy", 1, "--Jz", .5, "--neighbors", 2,
                  "--uniform-effective", "--g", .3, "--mass", 1.7,
                  "--omega-ph", .8, "--uniform-grid", 64,
                  "--extra-q", "zero,0,0", "--extra-q", "off,.137,.271")
        for geometry in ((6, 6, 0), (12, 3, 0), (12, 3, 3)):
            with self.subTest(geometry=geometry):
                report, cpp = self.generate_and_parse("g%d_%d_%d" % geometry, *geometry, *common)
                zero = next(point for point in report["q_file_measurement_points"]
                            if "user_zero" in point["aliases"])
                off = next(point for point in report["q_file_measurement_points"]
                           if "user_off" in point["aliases"])
                self.assertTrue(zero["on_grid"])
                self.assertFalse(off["on_grid"])
                self.assertTrue(report["uniform_background"]["converged"])
                self.assertEqual(len(cpp), report["q_file_deduplication"]["written_unique_points"])

    def test_duplicate_source_cases_are_cleanly_aliased(self):
        report, _ = self.generate_and_parse(
            "nearest_pm", 2, 2, 0, "--Jz", .05, "--neighbors", 1)
        self.assertTrue(any("bare_Qp_g3" in p["aliases"] and "bare_Qm_g3" in p["aliases"]
                            for p in report["q_file_measurement_points"]))

        report, _ = self.generate_and_parse(
            "bare_kek", 6, 6, 0, "--Jz", 1, "--neighbors", 1)
        self.assertTrue(any(any(a.startswith("bare_Q") for a in p["aliases"])
                            and any(a.startswith("kek_K") for a in p["aliases"])
                            for p in report["q_file_measurement_points"]))

        report, _ = self.generate_and_parse(
            "bare_uniform", 6, 6, 0, "--uniform-effective", "--g", 0,
            "--spring-k", 1, "--uniform-grid", 64, "--neighbors", 1)
        self.assertTrue(any(any(a.startswith("bare_Q") for a in p["aliases"])
                            and any(a.startswith("uniform_effective_Q") for a in p["aliases"])
                            for p in report["q_file_measurement_points"]))

        report, _ = self.generate_and_parse(
            "special_minus", 12, 3, 3, "--neighbors", 1,
            "--extra-q", "special_Q,.5,.5", "--extra-q", "special_minus,-.5,-.5",
            "--extra-q", "off_left,.01,.01", "--extra-q", "off_right,-.01,-.01")
        special = next(p for p in report["q_file_measurement_points"]
                       if "user_special_Q" in p["aliases"])
        self.assertIn("user_special_minus", special["aliases"])
        self.assertTrue(special["on_grid"])
        left = next(p for p in report["q_file_measurement_points"] if "user_off_left" in p["aliases"])
        right = next(p for p in report["q_file_measurement_points"] if "user_off_right" in p["aliases"])
        self.assertNotEqual(left["index"], right["index"])

    def test_cpp_parser_still_rejects_dirty_duplicate_input(self):
        with tempfile.TemporaryDirectory(prefix="vbs-q-dirty-") as directory:
            path = Path(directory) / "duplicate.txt"
            path.write_text("first .25 0\nsecond 1.25 0\n")
            parsed = subprocess.run([str(self.validator), "--inspect-q-file", "6", "6", "0", str(path)],
                                    text=True, capture_output=True)
            self.assertNotEqual(parsed.returncode, 0)
            self.assertIn("duplicate momentum", parsed.stderr)


if __name__ == "__main__":
    unittest.main()
