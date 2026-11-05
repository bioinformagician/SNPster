import copy
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pedigree_benchmark/scripts"))
from flare_nextflow_tasks import build_jobs, crossed_theta, SOURCES


class NextflowFlareTests(unittest.TestCase):
    def test_theta_transfer_uses_named_axes_and_preserves_other_parameters(self):
        names = list(SOURCES)
        donor_names = list(reversed(SOURCES))
        model = dict(ancestries=names, panels=names[1:] + names[:1], p=np.full((7, 7), 1/7),
                     theta=np.full((7, 7), 1e-5), rho=np.arange(1, 8), mu=np.full(7, 1/7), generation=10.)
        donor = dict(ancestries=donor_names, panels=donor_names[2:] + donor_names[:2],
                     theta=np.arange(49).reshape(7, 7) * 1e-6 + 1e-6)
        old = copy.deepcopy(model)
        result = crossed_theta(model, donor)
        for i, a in enumerate(model["ancestries"]):
            for j, b in enumerate(model["panels"]):
                self.assertEqual(result["theta"][i, j], donor["theta"][
                    donor["ancestries"].index(a), donor["panels"].index(b)])
        for key in ("p", "rho", "mu", "generation"):
            np.testing.assert_array_equal(result[key], old[key])
        np.testing.assert_array_equal(model["theta"], old["theta"])
        with self.assertRaises(ValueError):
            crossed_theta(model, dict(donor, ancestries=["bad"] + donor_names[1:]))

    def protocol(self):
        return dict(recovery=dict(candidates=[
            "finecap16_fixed", "finecap16_em", "size_control_fixed", "size_control_em"]),
                    stages=["inherited_phase", "array_phase_only"],
                    theta_controls=dict(panels=["lonni", "frederik"], references=["full", "finecap16", "size_control"]))

    def test_job_cartesian_coverage_no_reuse_or_evaluation(self):
        people = [dict(sample_id=f"train{i}", track="pure", source_group=SOURCES[i % 7]) for i in range(316)]
        jobs = build_jobs(dict(people=people), [], self.protocol())
        self.assertEqual(len(jobs), 6320)
        self.assertEqual(sum(j["kind"] == "recovery" for j in jobs), 2528)
        self.assertEqual(sum(j["kind"] == "theta" for j in jobs), 3792)
        self.assertEqual(len({j["id"] for j in jobs}), 6320)
        self.assertTrue(all(j["em"] for j in jobs if j["kind"] == "theta"))
        with self.assertRaises(ValueError):
            build_jobs(dict(people=people), [jobs[0]], self.protocol())

    def test_smoke_is_training_only_and_order_invariant(self):
        people = [dict(sample_id=f"SAS{i}", track="pure", source_group="SAS") for i in range(3)]
        people += [dict(sample_id=f"child{i}", track="pedigree", source_group="mixture") for i in range(3)]
        truth = {"22": {f"child{i}": ([1, 0, 0, 0, 0, 0, 0] if i < 2 else [.5, .5, 0, 0, 0, 0, 0])
                        for i in range(3)}}
        selection = dict(people=people, truth=truth)
        a = build_jobs(selection, [], self.protocol(), True)
        b = build_jobs(dict(selection, people=list(reversed(people))), [], self.protocol(), True)
        self.assertEqual(a, b)
        self.assertEqual(len(a), 40)
        self.assertIn("child2", {r["sample_id"] for r in a})


if __name__ == "__main__":
    unittest.main()
