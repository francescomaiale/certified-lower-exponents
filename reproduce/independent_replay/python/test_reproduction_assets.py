"""Check original GM assets and the floating-point Table S5 provenance/export."""
import hashlib
import json
from decimal import Decimal, ROUND_CEILING
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]


class ReproductionAssetsTests(unittest.TestCase):
    def test_original_gm_files_and_required_seeds(self):
        folder = ROOT / 'results/gm_original'
        manifest = json.loads((folder / 'MANIFEST.json').read_text())
        files = {r['path']: r for r in manifest['files']}
        self.assertEqual(len(files), 284)
        for name, record in files.items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((folder/name).read_bytes()).hexdigest(), record['sha256'])
        for k in (1, 3, 4, 6, 8, 9, 10):
            name = f'switching/sw10_{k}_0.1_E_V.mat'
            self.assertEqual(files[name]['shape'][0], 10)
            self.assertGreater(files[name]['shape'][1], 0)

    def test_table9_word_column_matches_original_logs(self):
        folder = ROOT / 'results/larger_samples'
        records = json.loads((folder/'scaling_word_estimates.json').read_text())['records']
        info = json.loads((ROOT/'data/scaling_pairs_journal_info.json').read_text())
        expected = {(r['id'], h) for r in info for h in ('0.1', '0.05')}
        self.assertEqual(len(records), 66)
        self.assertEqual({(r['pair'], r['h']) for r in records}, expected)
        tex = (ROOT/'results/larger_samples_tables.tex').read_text().split('% Table 9 (tab:scalingnew)')[1]
        rows = {}
        for line in tex.splitlines():
            if not re.match(r'\d+ & ', line):
                continue
            self.assertTrue(line.endswith('\\' * 2) and not line.endswith('\\' * 3))
            cells = [c.strip() for c in line.split('&')]
            self.assertEqual(len(cells), 9)
            kind, k = cells[1].split()
            rows[f'd{cells[0]}_{kind}{k}'] = re.search(r'[+-]?\d+\.\d{6}', cells[4])[0]
        self.assertEqual(set(rows), {r['id'] for r in info})
        excerpts = (folder/'scaling_word_log_excerpts.txt').read_text()
        for record in records:
            m = re.search(r'best explicit word \(k = (\d+)\): g\(h\) <= ([+\-0-9.]+)', record['source_line'])
            self.assertIsNotNone(m)
            self.assertEqual(int(m[1]), record['word_length'])
            self.assertEqual(m[2], record['g_up_logged'])
            self.assertIn(f"[{record['source_log']}]\n{record['source_line']}", excerpts)
            if record['h'] == '0.1':
                rounded = Decimal(record['g_up_logged']).quantize(Decimal('.000001'), rounding=ROUND_CEILING)
                self.assertEqual(Decimal(rows[record['pair']]), rounded)


if __name__ == '__main__':
    unittest.main()
