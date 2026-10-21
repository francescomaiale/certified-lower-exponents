"""Shared helpers of the Python replays: the merged data of ../../data/*.json, keyed by pair name.

The replay functions derive from the paper's Python scripts, with explicit input and certificate checks.
They share no code with the MATLAB certifiers.
"""
import atexit, json, os, tempfile
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "data")
MAP_FILES = ("test_pairs.json", "scaling_pairs.json", "switching_pairs.json", "bcv13.json",
             "scaling_pairs_journal.json", "switching_pairs_journal.json")   # the last two: larger samples, if present


def read_json(path, **kwargs):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh, **kwargs)


def check_family(As):
    """Reject data outside the nonempty, square Metzler setting of the paper."""
    if not As or not As[0]:
        raise ValueError("empty matrix family")
    d = len(As[0])
    for A in As:
        if len(A) != d or any(len(row) != d for row in A):
            raise ValueError("matrices must be square and have the same size")
        if any(A[r][c] < 0 for r in range(d) for c in range(d) if r != c):
            raise ValueError("matrices must be Metzler")


def vector(xs, d):
    if not isinstance(xs, list) or len(xs) != d:
        raise ValueError("wrong witness dimension")
    return [Fraction(x) for x in xs]


def word(value, m):
    """Legacy digit string, or a list of one-based integer mode indices."""
    if isinstance(value, str):
        if not value or any(c not in "123456789" for c in value):
            raise ValueError("invalid word")
        value = [int(c) for c in value]
    if not isinstance(value, list) or not value or any(type(i) is not int or not 1 <= i <= m for i in value):
        raise ValueError("invalid word")
    return [i - 1 for i in value]


def merged_cases_file():
    """Write the union of the pair maps to a temporary JSON file (numbers verbatim) and return its path."""
    parts = []
    for f in MAP_FILES:
        if not os.path.exists(os.path.join(DATA, f)):
            continue
        with open(os.path.join(DATA, f), encoding="utf-8") as fh:
            txt = fh.read().strip()
        assert txt[0] == "{" and txt[-1] == "}"
        parts.append(txt[1:-1].strip())
    fd, path = tempfile.mkstemp(suffix=".json", prefix="paperB_cases_")
    atexit.register(lambda: os.path.exists(path) and os.remove(path))
    with os.fdopen(fd, "w") as fh:
        fh.write("{" + ",\n".join(parts) + "}\n")
    read_json(path)                            # well formed
    return path
