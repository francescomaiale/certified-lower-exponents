"""Numerical check of the rate-constant theorem (easy half) for paper B, Appendix C / Table 5.

Claim (two modes, interior chattering minimizer alpha*, S = S_alpha* irreducible, u^T S = eta u^T, S v = eta v,
u^T v = 1, F = A_1 - A_2, Q = alpha*_1 alpha*_2 u^T F^2 v): along a balanced word w of period q with frequency beta,
T = q h,
    log rho(B_w) / T = eta_beta - (h/2) u_b^T (sum_i beta_i A_i^2) v_b + O(hT + h^2)      (B_i = I + h A_i)
    log rho(R_w) / T = eta_beta + (h/2) u_b^T (sum_i beta_i A_i^2) v_b + O(hT + h^2)      (R_i = (I - h A_i)^-1)
so the word bounds g_w = (rho(B_w)^(1/q) - 1)/h >= g(h) and G_w = (1 - rho(R_w)^(-1/q))/h >= G(h) satisfy
    g_w = c - (Q/2) h + O(h^(4/3)),   G_w = c + (Q/2) h + O(h^(4/3))   for q = ceil(h^(-2/3)).
Usage: python3 check_rate_theorem.py   (writes ../results/check_rate_theorem.md; about 2 minutes)
"""
import json, math, os
import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def load(name):
    return {k: [np.array(M, float) for M in v] for k, v in json.load(open(os.path.join(DATA, name))).items()}


def eta(M):
    return max(np.linalg.eigvals(M).real)


def perron(S):
    w, V = np.linalg.eig(S); k = np.argmax(w.real); v = np.real(V[:, k]); v = v * np.sign(v.sum())
    w2, U = np.linalg.eig(S.T); k2 = np.argmax(w2.real); u = np.real(U[:, k2]); u = u / (u @ v)
    return w[k].real, u, v


def chatter(A1, A2):
    r = minimize_scalar(lambda a: eta(a * A1 + (1 - a) * A2), bounds=(0, 1), method="bounded",
                        options={"xatol": 1e-13})
    a = r.x
    return a, eta(a * A1 + (1 - a) * A2)


def constants(A1, A2, a):
    S = a * A1 + (1 - a) * A2; e, u, v = perron(S); F = A1 - A2
    QB = a * (1 - a) * (u @ F @ F @ v)
    Q_alt = u @ (a * A1 @ A1 + (1 - a) * A2 @ A2) @ v - e ** 2          # equivalent form in Theorem 3.12
    QA = 0.5 * (u @ (F @ (F @ S - S @ F) - (F @ S - S @ F) @ F) @ v)   # paper A's high-frequency quantity
    h2 = 1e-4                                                           # eta'' by central differences
    d2 = (eta((a + h2) * A1 + (1 - a - h2) * A2) - 2 * eta(S) + eta((a - h2) * A1 + (1 - a + h2) * A2)) / h2 ** 2
    return QB, Q_alt, QA, d2


def word(q, p):
    """lower mechanical (Christoffel) word of length q with p letters 1: discrepancy < 1."""
    return [1 if (j * p) // q - ((j - 1) * p) // q == 1 else 2 for j in range(1, q + 1)]


def rho(M):
    return max(abs(np.linalg.eigvals(M)))


def word_bounds(A1, A2, h, w):
    d = A1.shape[0]; I = np.eye(d)
    B = {1: I + h * A1, 2: I + h * A2}; R = {1: np.linalg.inv(I - h * A1), 2: np.linalg.inv(I - h * A2)}
    PB = I.copy(); PR = I.copy()
    for x in w:
        PB = B[x] @ PB; PR = R[x] @ PR
    q = len(w); rB, rR = rho(PB), rho(PR)
    return (rB ** (1 / q) - 1) / h, (1 - rR ** (-1 / q)) / h, math.log(rB) / (q * h), math.log(rR) / (q * h)


def main():
    out = []; say = lambda s="": (print(s), out.append(s))
    test = load("test_pairs.json"); sc = load("scaling_pairs.json")
    pairs = [("planar", test["planar"]), ("case 1", test["case1"])] + [(k, sc[k]) for k in
             ["d5_dense1", "d5_dense2", "d5_dense3", "d5_sparse1", "d5_sparse2", "d5_sparse3",
              "d10_dense1", "d10_dense2", "d10_dense3", "d10_sparse1", "d10_sparse2", "d10_sparse3"]]
    table11 = [1.4780, 0.5480, 0.1785, 0.5888, 0.1975, 1.3307, 0.2024, 2.3887, 0.0447, 0.3178, 0.0770, 0.7134,
               0.4656, 1.5941]
    say("# Rate-constant theorem (easy half): numerical check")
    say()
    say("## 1. Q against Table 5")
    say()
    say("| pair | alpha*_1 | c | Q = a(1-a) u'F^2v | Theorem 3.12 form | Table 5 | eta'' | paper A's Q |")
    say("|---|---|---|---|---|---|---|---|")
    info = {}
    for (name, (A1, A2)), q11 in zip(pairs, table11):
        a, c = chatter(A1, A2); QB, Qalt, QA, d2 = constants(A1, A2, a); info[name] = (a, c, QB)
        say(f"| {name} | {a:.4f} | {c:+.6f} | {QB:.4f} | {Qalt:.4f} | {q11:.4f} | {d2:.3g} | {QA:+.4g} |")
    say()
    say("## 2. Balanced-word bounds, q = ceil(h^(-2/3)), p = round(q alpha*_1)")
    say()
    say("Columns: (c - g_w)/h and (G_w - c)/h, which the theorem says tend to Q/2; their sum tends to Q; the")
    say("deviations from Q/2 divided by h^(1/3) stay bounded.")
    hs = [0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125]
    for name in ["planar", "case 1", "d5_dense2", "d5_sparse3", "d10_sparse1"]:
        A1, A2 = dict(pairs)[name]; a, c, QB = info[name]
        say(); say(f"**{name}** (Q/2 = {QB / 2:.5f})"); say()
        say("| h | q | p/q | (c - g_w)/h | (G_w - c)/h | sum | dev. expl./h^(1/3) | dev. impl./h^(1/3) |")
        say("|---|---|---|---|---|---|---|---|")
        for h in hs:
            q = math.ceil(h ** (-2 / 3)); p = min(q - 1, max(1, round(q * a)))
            gw, Gw, _, _ = word_bounds(A1, A2, h, word(q, p))
            e1 = (c - gw) / h; e2 = (Gw - c) / h
            say(f"| {h:g} | {q} | {p}/{q} | {e1:.5f} | {e2:.5f} | {e1 + e2:.5f} | {(e1 - QB / 2) / h ** (1 / 3):+.3f} |"
                f" {(e2 - QB / 2) / h ** (1 / 3):+.3f} |")
    say()
    say("## 2b. Balanced-word bounds with a fixed period, q = floor(T1/h), T1 = 0.5 (the O(h^2) version)")
    say()
    say("Deviations of (c - g_w)/h and (G_w - c)/h from Q/2, divided by h: bounded if the bounds are c -+ (Q/2)h + O(h^2).")
    for name in ["planar", "case 1", "d5_dense2", "d5_sparse3", "d10_sparse1", "d10_sparse3"]:
        A1, A2 = dict(pairs)[name]; a, c, QB = info[name]
        say(); say(f"**{name}** (Q/2 = {QB / 2:.5f})"); say()
        say("| h | q | p/q | (c - g_w)/h | (G_w - c)/h | dev. expl./h | dev. impl./h |")
        say("|---|---|---|---|---|---|---|")
        for h in hs:
            q = int(0.5 / h); p = min(q - 1, max(1, round(q * a)))
            gw, Gw, _, _ = word_bounds(A1, A2, h, word(q, p))
            e1 = (c - gw) / h; e2 = (Gw - c) / h
            say(f"| {h:g} | {q} | {p}/{q} | {e1:.5f} | {e2:.5f} | {(e1 - QB / 2) / h:+.3f} | {(e2 - QB / 2) / h:+.3f} |")
    say()
    say("## 3. The lemma at a fixed frequency: remainder / h^2")
    say()
    say("For beta = p/q fixed, log rho(B_w)/T - [eta_beta - (h/2) u_b'(sum beta_i A_i^2) v_b] and the implicit analogue,")
    say("divided by h^2: bounded, and independent of q, if the lemma holds with remainder O(h^2).")
    say()
    say("| pair | q | h | T | explicit | implicit |")
    say("|---|---|---|---|---|---|")
    for name in ["planar", "case 1", "d10_sparse3"]:
        A1, A2 = dict(pairs)[name]; a = info[name][0]
        for q in (8, 32, 128):
            p = max(1, min(q - 1, round(q * a))); b = p / q; S = b * A1 + (1 - b) * A2; e, u, v = perron(S)
            m2 = u @ (b * A1 @ A1 + (1 - b) * A2 @ A2) @ v
            for h in (0.004, 0.002, 0.001):
                _, _, lB, lR = word_bounds(A1, A2, h, word(q, p)); T = q * h
                say(f"| {name} | {q} | {h:g} | {T:g} | {(lB - (e - h / 2 * m2)) / (h * h):+.3f} |"
                    f" {(lR - (e + h / 2 * m2)) / (h * h):+.3f} |")
    say()
    say("## 4. Q as the mean of paper A's Floquet function h(w) = Re u'F(eta + iw - S)^-1 F v")
    say()
    say("Identity to check: Q = (2 a(1-a)/pi) * int_0^inf h(w) dw.")
    say()
    say("| pair | Q | (2a(1-a)/pi) int h | min_w h(w) |")
    say("|---|---|---|---|")
    from scipy.integrate import quad
    for name in ["planar", "case 1", "d5_sparse1", "d10_dense2", "d10_sparse3"]:
        A1, A2 = dict(pairs)[name]; a, c, QB = info[name]
        S = a * A1 + (1 - a) * A2; e, u, v = perron(S); F = A1 - A2; d = S.shape[0]
        hw = lambda w: np.real(u @ F @ np.linalg.solve((e + 1j * w) * np.eye(d) - S, F @ v))
        val = sum(quad(hw, lo, hi, limit=400)[0] for lo, hi in [(0, 1), (1, 10), (10, 100), (100, 1e4)]) + \
              0.5 * (u @ (F @ (F @ S - S @ F) - (F @ S - S @ F) @ F) @ v) / 1e4     # tail ~ Q_A / w^2
        ws = np.concatenate([np.logspace(-3, 0, 300), np.linspace(1, 5, 401), np.logspace(0.7, 4, 400)])
        say(f"| {name} | {QB:.6f} | {2 * a * (1 - a) / math.pi * val:.6f} | {min(hw(w) for w in ws):+.4g} |")
    say()
    say("## 5. Sign of Q on random pairs with an interior chattering minimizer")
    say()
    rng = np.random.default_rng(1)
    for d, kind in [(3, "sparse"), (3, "dense"), (5, "sparse"), (5, "dense")]:
        n = nneg = nA = nboth = 0
        while n < 400:
            def draw():
                if kind == "dense":
                    M = rng.uniform(0, 2 / (d - 1), (d, d))
                else:
                    M = (rng.random((d, d)) < 0.5) * rng.lognormal(0, 1.5, (d, d)) * 2.46 / (0.5 * (d - 1) * math.exp(1.125))
                np.fill_diagonal(M, rng.uniform(-3, 1, d)); return M
            A1, A2 = draw(), draw()
            a, c = chatter(A1, A2)
            if not 1e-3 < a < 1 - 1e-3:
                continue
            S = a * A1 + (1 - a) * A2; ev = np.sort(np.linalg.eigvals(S).real)
            if ev[-1] - ev[-2] < 1e-6:
                continue
            n += 1; QB, _, QA, _ = constants(A1, A2, a)
            nneg += QB < 0; nA += QA < 0; nboth += (QB < 0) and (QA < 0)
        say(f"- d = {d}, {kind}: {n} pairs with an interior minimizer; Q < 0 in {nneg}, paper A's high-frequency test < 0 in {nA}, both in {nboth}")
    open(os.path.join(HERE, "..", "results", "check_rate_theorem.md"), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
