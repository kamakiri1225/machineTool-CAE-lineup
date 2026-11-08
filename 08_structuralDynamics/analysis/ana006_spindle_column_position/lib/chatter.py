"""相対 FRF（モード和）とフライス安定限界（ゼロ次近似、ana005 5.2 と同じ式）"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import params as P

B0 = P.TEETH / 4 * np.array([[P.KR, 1.0], [-1.0, P.KR]])
GRID = np.linspace(5, 3000, 24000)
RPM = np.linspace(500, 16000, 600)


def damping(fr):
    wr = 2 * np.pi * fr
    return wr, P.ALPHA / (2 * wr) + P.BETA * wr / 2


def G_of(f, fr, d):
    """G(f) = Σ d_r d_r^T / (ω_r² − ω² + 2iζ_r ω_r ω)   d: (2 or 3, nm)"""
    wr, zr = damping(fr); w = 2 * np.pi * np.atleast_1d(f)
    den = wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]
    return np.einsum("ir,jr,wr->wij", d, d, 1 / den)


def boundary(interp, rpm, grid, Kc=P.KT, teeth=P.TEETH):
    res, crit = [], []
    for n in rpm:
        tau = 60 / (teeth * n); h = interp(grid) * (1 - np.exp(-2j * np.pi * grid * tau))
        idx = np.where(h.imag[:-1] * h.imag[1:] < 0)[0]; cand = []
        for i in idx:
            f = brentq(lambda q: (interp(q) * (1 - np.exp(-2j * np.pi * q * tau))).imag, grid[i], grid[i + 1], xtol=1e-9)
            v = interp(f) * (1 - np.exp(-2j * np.pi * f * tau))
            if v.real < -1e-12:
                b = -1 / (Kc * v.real)
                if b < 1e6:
                    cand.append((b, f))
        b, f = min(cand) if cand else (np.nan, np.nan); res.append(b); crit.append(f)
    return np.array(res), np.array(crit)


def lobes(fr, d):
    Gg = G_of(GRID, fr, d[:2]); lam = np.linalg.eigvals(Gg @ B0)
    for i in range(1, len(lam)):
        if np.sum(abs(lam[i] - lam[i - 1])) > np.sum(abs(lam[i, ::-1] - lam[i - 1])):
            lam[i] = lam[i, ::-1]
    br = [boundary(interp1d(GRID, lam[:, j], bounds_error=True), RPM, GRID) for j in (0, 1)]
    v = np.stack([b[0] for b in br], 1); fq = np.stack([b[1] for b in br], 1)
    w = np.argmin(np.where(np.isfinite(v), v, np.inf), 1)
    return v[np.arange(len(RPM)), w], fq[np.arange(len(RPM)), w], Gg


def mode_share(fc, fr, d):
    """びびり周波数 fc での揺れ Q を各モードに分けた割合（ana005 5.3 と同じ）"""
    wr, zr = damping(fr); wc = 2 * np.pi * fc; den = wr ** 2 - wc ** 2 + 2j * zr * wr * wc
    G = np.einsum("ir,jr,r->ij", d[:2], d[:2], 1 / den); lam, V = np.linalg.eig(G @ B0)
    v = V[:, np.argmax((-1 / lam).real)]; F = B0 @ v; eta = (d[:2].T @ F) / den
    s = np.abs(eta) * np.linalg.norm(d[:2], axis=0); return s / s.sum()
