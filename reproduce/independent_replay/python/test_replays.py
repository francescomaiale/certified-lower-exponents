"""Regression tests for the proof prerequisites and malformed-certificate rejection.

Run: python3 -m unittest discover -s reproduce/independent_replay/python -p 'test_*.py'
The scalar counterexamples have known exponents, so rejection is mathematically necessary.
"""
import contextlib
import copy
from decimal import Decimal, localcontext
from fractions import Fraction as Fr
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import common
import replay_upper
import replay_best_upper
import replay_lower
import replay_cones
import replay_faces
import run_all


class ReplayTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.data = self.root / 'data.json'
        self.cert = self.root / 'certificate.json'

    def upper(self, C, raw):
        self.data.write_text(json.dumps({'scalar': raw}))
        self.cert.write_text(json.dumps(C))
        with patch.object(replay_upper, 'CASES', {'test': str(self.data)}):
            return replay_upper.replay_one(str(self.cert))[0]

    def best(self, C, raw):
        self.data.write_text(json.dumps({'scalar': raw}))
        self.cert.write_text(json.dumps(C))
        with patch.object(replay_best_upper, 'CASES', str(self.data)), contextlib.redirect_stdout(io.StringIO()):
            return replay_best_upper.replay_file(str(self.cert))

    @staticmethod
    def upper_cert():
        rd = {'z': [['1'], ['1']], 'chain': [['1'], ['1']], 'r': '1'}
        return {'case': 'scalar', 'source': 'test', 'h': '0.1', 'word': '1', 'k': 1,
                'printed': '0', 'printed_exact': '0',
                'versions': {v: copy.deepcopy(rd) for v in ('decimal', 'binary64')}}

    @staticmethod
    def best_cert(schedule=False):
        rd = ({'modes': [0], 'durations': ['1'], 'x_sched': ['1']} if schedule
              else {'alpha': '1', 'x': ['1']})
        U = '1/1000' if schedule else '0'
        return {'id': 'scalar', 'U_c': [] if schedule else U, 'U_s': U if schedule else [],
                'U_best': U, 'printed': U,
                'readings': {v: copy.deepcopy(rd) for v in ('decimal', 'binary64')}}

    def test_valid_upper_and_single_mode(self):
        C = self.upper_cert()
        self.assertTrue(self.upper(C, [[[0]], [[0]]]))
        for rd in C['versions'].values():
            rd['z'] = [['1']]
        self.assertTrue(self.upper(C, [[[0]]]))

    def test_non_mmatrix_counterexample(self):
        # A_1 = A_2 = 10, lambda = 10; the old replay accepted the false bound lambda <= 0.
        C = self.upper_cert()
        C.update(h='0.2', word='11', k=2)
        for rd in C['versions'].values():
            rd.update(z=[], chain=[['1'], ['-1'], ['1']])
        self.assertFalse(self.upper(C, [[[10]], [[10]]]))

    def test_upper_prerequisites(self):
        for field, value in [('h', '0'), ('h', '-0.1'), ('word', ''), ('word', '0'),
                             ('word', '3'), ('word', [1.5]), ('k', 2)]:
            with self.subTest(field=field, value=value):
                C = self.upper_cert(); C[field] = value
                self.assertFalse(self.upper(C, [[[0]], [[0]]]))
        for field, value in [('z', []), ('z', [['1']]), ('z', [['1', '1'], ['1']]),
                             ('chain', []), ('chain', [['0'], ['1']]), ('r', '-1'), ('r', '0')]:
            with self.subTest(field=field, value=value):
                C = self.upper_cert(); C['versions']['decimal'][field] = value
                self.assertFalse(self.upper(C, [[[0]], [[0]]]))

    def test_numeric_word_supports_ten_modes(self):
        C = self.upper_cert(); C['word'] = [10]
        for rd in C['versions'].values():
            rd['z'] = [['1'] for _ in range(10)]
        self.assertTrue(self.upper(C, [[[0]] for _ in range(10)]))

    def test_mixture_counterexample(self):
        C = self.best_cert()
        self.assertTrue(self.best(C, [[[0]], [[1]]]))
        # alpha=2 produces -1, although every actual switching law has exponent >= 0.
        for rd in C['readings'].values():
            rd['alpha'] = '2'
        C.update(U_c='-1', U_best='-1', printed='-1')
        self.assertFalse(self.best(C, [[[0]], [[1]]]))

    def test_schedule_counterexample(self):
        C = self.best_cert(schedule=True)
        self.assertTrue(self.best(C, [[[0]], [[1]]]))
        # Old zip() used exp(A*1) but divided by T=2: a false 0.75 bound for lambda=1.
        for rd in C['readings'].values():
            rd['durations'] = ['1', '1']
        C.update(U_s='3/4', U_best='3/4', printed='0.75')
        self.assertFalse(self.best(C, [[[1]], [[1]]]))

    def test_schedule_prerequisites(self):
        for field, value in [('modes', []), ('modes', [-1]), ('modes', [0.5]), ('modes', [2]),
                             ('durations', ['0']), ('durations', ['-1']), ('x_sched', ['0']), ('x_sched', [])]:
            with self.subTest(field=field, value=value):
                C = self.best_cert(schedule=True); C['readings']['decimal'][field] = value
                self.assertFalse(self.best(C, [[[0]], [[1]]]))

    def test_printed_best_bound(self):
        C = self.best_cert(); C['printed'] = '-1'
        self.assertFalse(self.best(C, [[[0]], [[1]]]))

    def test_large_negative_exponential(self):
        with localcontext() as ctx:
            ctx.prec = 180
            for y in [-1, -81, -82, -100, -200]:
                with self.subTest(y=y):
                    q = replay_best_upper.exp_lower(Fr(y))
                    self.assertGreater(q, 0)
                    self.assertLessEqual(Decimal(q.numerator) / Decimal(q.denominator), Decimal(y).exp())

    def test_empty_and_unknown_selection(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(run_all.main([str(self.root / 'missing')]))
            self.assertFalse(run_all.main([str(self.root)]))
            self.cert.write_text('{}')
            self.assertFalse(run_all.main([str(self.root)]))
            for mod in (replay_lower, replay_upper, replay_best_upper, replay_cones, replay_faces):
                with self.subTest(module=mod.__name__):
                    self.assertFalse(mod.main([]))

    def test_existing_lower_negative_step_and_orientation(self):
        C = json.loads((Path(common.ROOT) / 'tests/fixtures/python_certificates/lower_case148_0.1_rowsumN_trivial.json').read_text())
        for field, value in [('h', '-0.1'), ('orientation', 'invalid')]:
            with self.subTest(field=field):
                T = copy.deepcopy(C); T[field] = value
                self.cert.write_text(json.dumps(T))
                self.assertFalse(replay_lower.replay_one(str(self.cert))[0])

    def test_cone_support_classes_and_empty_faces(self):
        paths = list((Path(common.ROOT) / 'certificates/comp38').glob('*.json'))
        cone = next(json.loads(p.read_text()) for p in paths if 'h1' in json.loads(p.read_text()))
        cone['classes'] = []
        self.cert.write_text(json.dumps(cone))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(replay_cones.replay_file(str(self.cert)))
            self.cert.write_text(json.dumps({'case': cone['case'], 'orientations': []}))
            self.assertFalse(replay_faces.replay_file(str(self.cert)))


if __name__ == '__main__':
    unittest.main()
