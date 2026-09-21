#!/usr/bin/env python3
"""Load x gain hazard on a subsystem graph, and what a demographic schedule can see.

Local damage is an age-linear ODE. The hazard map is either additive in the
mean load or exponential in it (a scalar load x gain). A second gain sits in
the damage rates themselves, as coupling to the other nodes. Fisher ranks and
profile likelihoods ask which of those parameters a synthetic age-specific
hazard schedule can separate.

The schedule is not a download from the Human Mortality Database. Exposure is
a declared synthetic curve. Node names are indices, not organs and not diagnoses.

Research computation only. Not a medical device, a dose, or a claim of
age reversal.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parent
FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)

SEED = 20260921

# Cuts declared with the design, before the scientific reading of the output.
# Eigenvalues below ABS_EIG_FLOOR are the floating-point null cluster.
# A relative cut on the raw matrix is not the rank: alpha and a per-year
# slope do not share units, and that cut mis-labels a real direction.
ABS_EIG_FLOOR = 1e-4
COLUMN_SV_TOL = 1e-8
CHI2_95_1 = 3.841
FLAT_CHI2 = 0.5
NEAR_LINEAR_R2 = 0.995
NEAR_LINEAR_MAX_REL = 0.02

N = 6
A0 = 30.0
A1 = 95.0
AGES = np.arange(A0, A1 + 1.0, 1.0)
T = AGES - A0

# Generating damage rates and initial loads. Model units. Not a life table.
V = np.array([0.010, 0.012, 0.009, 0.011, 0.013, 0.008], dtype=float)
D0 = np.array([0.20, 0.18, 0.22, 0.19, 0.21, 0.17], dtype=float)
GAMMA = 8.0
MU_B = 2.0e-4
LOG_MU_B = float(np.log(MU_B))
VBAR = float(np.mean(V))
L0 = float(np.mean(D0))

# Synthetic exposure: shrinking risk set, not a national population.
E0 = 2.0e5
EXPOSURE_DECAY = 0.035
EXPOSURE_FLOOR = 500.0

SIGMA_L = 0.02
SIGMA_D = 0.02

# Damage-rate coupling values fixed before the run.
G_SCAN = np.array([0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08], dtype=float)
G_PROFILE_TRUTHS = (0.0, 0.02)

GAMMA_FACTORS = np.geomspace(0.25, 4.0, 17)
GAMMA_FINE_FACTORS = np.linspace(0.92, 1.08, 33)
BETA_FACTORS = np.linspace(0.97, 1.03, 31)
G_FINE = {
    0.0: np.linspace(-0.003, 0.003, 31),
    0.02: np.linspace(0.016, 0.024, 33),
}


def exposure(t: np.ndarray) -> np.ndarray:
    return np.maximum(E0 * np.exp(-EXPOSURE_DECAY * t), EXPOSURE_FLOOR)


def sigma_log_from_mu(mu: np.ndarray, expo: np.ndarray) -> np.ndarray:
    deaths = np.maximum(mu * expo, 1e-8)
    return 1.0 / np.sqrt(deaths)


def damage_linear(t: np.ndarray) -> np.ndarray:
    return D0[None, :] + np.outer(t, V)


def coupling_matrix() -> np.ndarray:
    """Each node feels the mean of the others. Eigenvalue 1 on the mean mode."""
    return (np.ones((N, N)) - np.eye(N)) / (N - 1)


def mean_damage_scalar(t: np.ndarray, v: float, load0: float, g: float) -> np.ndarray:
    """Solution of dm/da = v + g m, m(a0) = load0. Stable at g = 0."""
    t = np.asarray(t, dtype=float)
    gt = g * t
    if abs(g) < 1e-10:
        series = t * (1.0 + gt / 2.0 + gt**2 / 6.0 + gt**3 / 24.0 + gt**4 / 120.0)
        return np.exp(gt) * load0 + v * series
    return np.exp(gt) * load0 + v * (np.expm1(gt) / g)


def damage_network(t: np.ndarray, g: float) -> np.ndarray:
    if abs(g) < 1e-12:
        return damage_linear(t)
    mate = g * coupling_matrix()
    mate_inv = np.linalg.inv(mate)
    eye = np.eye(N)
    out = np.empty((t.size, N), dtype=float)
    for i, ti in enumerate(t):
        big = expm(mate * float(ti))
        out[i] = big @ D0 + mate_inv @ ((big - eye) @ V)
    return out


def log_mu_gain(t: np.ndarray, log_mu_b: float, gamma: float, v: float, load0: float) -> np.ndarray:
    return log_mu_b + gamma * (load0 + v * t)


def mu_from_log(log_mu: np.ndarray) -> np.ndarray:
    return np.exp(log_mu)


def r2(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 0.0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def linearity(x: np.ndarray, y: np.ndarray) -> dict:
    design = np.column_stack([np.ones_like(x), x])
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    yhat = design @ coef
    increment = abs(float(coef[1]) * (float(x[-1]) - float(x[0])))
    max_abs = float(np.max(np.abs(y - yhat)))
    return {
        "r2": r2(y, yhat),
        "max_abs_resid": max_abs,
        "increment": increment,
        "max_rel_resid": max_abs / max(increment, 1e-12),
        "slope": float(coef[1]),
    }


def wls_poly(t: np.ndarray, y: np.ndarray, sigma: np.ndarray, deg: int) -> tuple[np.ndarray, float]:
    scale = 1.0 / sigma
    cols = [(t**k) * scale for k in range(deg + 1)]
    design = np.column_stack(cols)
    coef, _, _, _ = np.linalg.lstsq(design, y * scale, rcond=None)
    yhat = sum(float(coef[k]) * t**k for k in range(deg + 1))
    rss = float(np.sum(((y - yhat) / sigma) ** 2))
    return coef, rss


def gompertz_call(log_mu: np.ndarray, sigma: np.ndarray) -> dict:
    coef1, rss1 = wls_poly(T, log_mu, sigma, 1)
    coef2, rss2 = wls_poly(T, log_mu, sigma, 2)
    df = int(T.size - 2)
    crit = float(chi2.ppf(0.95, df))
    return {
        "alpha": float(coef1[0]),
        "beta": float(coef1[1]),
        "quadratic_c": float(coef2[2]),
        "rss_linear": rss1,
        "rss_quadratic": rss2,
        "lr_vs_quadratic": rss1 - rss2,
        "chi2_crit_df": crit,
        "df": df,
        "gompertz_like": bool(rss1 <= crit),
        "curvature_detectable": bool((rss1 - rss2) > CHI2_95_1),
    }


def node_linearity(damage: np.ndarray) -> dict:
    stats = [linearity(AGES, damage[:, i]) for i in range(damage.shape[1])]
    return {
        "min_r2": float(min(s["r2"] for s in stats)),
        "max_rel_resid": float(max(s["max_rel_resid"] for s in stats)),
        "near_linear": bool(
            min(s["r2"] for s in stats) >= NEAR_LINEAR_R2
            and max(s["max_rel_resid"] for s in stats) <= NEAR_LINEAR_MAX_REL
        ),
        "per_node_r2": [s["r2"] for s in stats],
    }


def ring_adjacency(c: float) -> np.ndarray:
    adj = np.zeros((N, N))
    for i in range(N):
        adj[i, (i + 1) % N] = c
    return adj


def star_adjacency(c: float) -> np.ndarray:
    adj = np.zeros((N, N))
    for i in range(1, N):
        adj[i, 0] = c
        adj[0, i] = c
    return adj


def felt(damage: np.ndarray, adj: np.ndarray) -> np.ndarray:
    return damage + damage @ adj.T


def graph_summary(name: str, adj: np.ndarray, compensate: bool) -> dict:
    damage = damage_linear(T)
    felt_load = felt(damage, adj)
    mean_felt = felt_load.mean(axis=1)
    # mean_felt = S0 + Sv * t, because damage is affine.
    design = np.column_stack([np.ones_like(T), T])
    coef, _, _, _ = np.linalg.lstsq(design, mean_felt, rcond=None)
    s0 = float(coef[0])
    sv = float(coef[1])
    beta = GAMMA * VBAR
    alpha = LOG_MU_B + GAMMA * L0
    if compensate:
        gamma = beta / sv
        log_mu_b = alpha - gamma * s0
    else:
        gamma = GAMMA
        log_mu_b = LOG_MU_B
    log_mu = log_mu_b + gamma * mean_felt
    target = log_mu_gain(T, LOG_MU_B, GAMMA, VBAR, L0)
    return {
        "name": name,
        "compensate": compensate,
        "gamma": float(gamma),
        "log_mu_b": float(log_mu_b),
        "mu_b": float(np.exp(log_mu_b)),
        "load_intercept": s0,
        "load_slope": sv,
        "max_abs_log_mu_gap": float(np.max(np.abs(log_mu - target))),
        "log_mu": log_mu,
        "felt": felt_load,
    }


def eig_report(fisher: np.ndarray) -> dict:
    fisher = 0.5 * (fisher + fisher.T)
    evals = np.linalg.eigvalsh(fisher)
    evals = np.array(evals[::-1], dtype=float)
    evals_clip = np.maximum(evals, 0.0)
    numerical = int(np.sum(evals_clip >= ABS_EIG_FLOOR))
    return {
        "eigenvalues_desc": evals_clip.tolist(),
        "min_unclipped": float(evals.min()) if evals.size else 0.0,
        "numerical_rank": numerical,
        "n_param": int(fisher.shape[0]),
        "condition_nonzero": (
            float(evals_clip[0] / evals_clip[numerical - 1])
            if numerical >= 1 and evals_clip[numerical - 1] > 0
            else None
        ),
    }


def column_spectrum(jac: np.ndarray, sigma: np.ndarray) -> dict:
    """Singular values after each weighted column is scaled to unit length."""
    weighted = jac / sigma[:, None]
    norms = np.linalg.norm(weighted, axis=0)
    norm_max = float(norms.max()) if norms.size else 0.0
    dead = norms <= COLUMN_SV_TOL * max(norm_max, 1.0)
    active = ~dead
    if not np.any(active):
        sv = np.array([], dtype=float)
    else:
        balanced = weighted[:, active] / norms[active]
        sv = np.linalg.svd(balanced, compute_uv=False)
    rank = int(np.sum(sv >= COLUMN_SV_TOL * sv[0])) if sv.size else 0
    rank += 0  # dead columns are excluded; they are not identifiable directions
    return {
        "column_rank": int(rank),
        "dead_columns": int(np.sum(dead)),
        "singular_values": sv.tolist(),
        "column_norms": norms.tolist(),
    }


def fisher_outer(jac: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    weighted = jac / sigma[:, None]
    return weighted.T @ weighted


def null_alignment(fisher: np.ndarray, vectors: list[np.ndarray]) -> list[float]:
    evals, evecs = np.linalg.eigh(0.5 * (fisher + fisher.T))
    mask = evals < ABS_EIG_FLOOR
    basis = evecs[:, mask]
    scores = []
    for vec in vectors:
        unit = vec / np.linalg.norm(vec)
        if basis.size == 0:
            scores.append(0.0)
        else:
            projected = basis @ (basis.T @ unit)
            scores.append(float(np.linalg.norm(projected)))
    return scores


def jac_mechanistic_4(gamma: float, v: float, load0: float, with_load: bool) -> tuple[np.ndarray, np.ndarray]:
    """Columns: log_mu_b, gamma, v, L0. Rows: log hazard, then optional mean load."""
    load = load0 + v * T
    j_h = np.column_stack(
        [
            np.ones_like(T),
            load,
            gamma * T,
            np.full_like(T, gamma),
        ]
    )
    sig_h = sigma_log_from_mu(mu_from_log(log_mu_gain(T, LOG_MU_B, GAMMA, VBAR, L0)), exposure(T))
    if not with_load:
        return j_h, sig_h
    j_l = np.column_stack(
        [
            np.zeros_like(T),
            np.zeros_like(T),
            T,
            np.ones_like(T),
        ]
    )
    return np.vstack([j_h, j_l]), np.concatenate([sig_h, np.full(T.size, SIGMA_L)])


def jac_nodes(with_load: str) -> tuple[np.ndarray, np.ndarray]:
    """14 columns: log_mu_b, gamma, v_1..v_6, D0_1..D0_6.

    with_load is 'hazard', 'mean', or 'nodes'.
    """
    load = L0 + VBAR * T
    n_age = T.size
    # hazard block
    j_h = np.zeros((n_age, 14))
    j_h[:, 0] = 1.0
    j_h[:, 1] = load
    for i in range(N):
        j_h[:, 2 + i] = GAMMA * T / N
        j_h[:, 8 + i] = GAMMA / N
    sig_h = sigma_log_from_mu(mu_from_log(log_mu_gain(T, LOG_MU_B, GAMMA, VBAR, L0)), exposure(T))
    blocks = [j_h]
    sigs = [sig_h]
    if with_load == "mean":
        j_m = np.zeros((n_age, 14))
        for i in range(N):
            j_m[:, 2 + i] = T / N
            j_m[:, 8 + i] = 1.0 / N
        blocks.append(j_m)
        sigs.append(np.full(n_age, SIGMA_L))
    elif with_load == "nodes":
        # age-major: for each age, N node rows. Column order matches parameters.
        j_n = np.zeros((n_age * N, 14))
        for k in range(n_age):
            for i in range(N):
                row = k * N + i
                j_n[row, 2 + i] = T[k]
                j_n[row, 8 + i] = 1.0
        blocks.append(j_n)
        sigs.append(np.full(n_age * N, SIGMA_D))
    return np.vstack(blocks), np.concatenate(sigs)


def cr_diag(fisher: np.ndarray) -> list[float] | None:
    report = eig_report(fisher)
    if report["numerical_rank"] < fisher.shape[0]:
        return None
    try:
        inv = np.linalg.inv(0.5 * (fisher + fisher.T))
    except np.linalg.LinAlgError:
        return None
    return np.maximum(np.diag(inv), 0.0).tolist()


def profile_gamma(
    y: np.ndarray,
    sigma: np.ndarray,
    z: np.ndarray | None,
    factors: np.ndarray | None = None,
) -> list[dict]:
    factors = GAMMA_FACTORS if factors is None else np.asarray(factors, dtype=float)
    def residuals(x: np.ndarray, gamma: float) -> np.ndarray:
        log_mu_b, v, load0 = (float(v_) for v_ in x)
        pred = log_mu_gain(T, log_mu_b, gamma, v, load0)
        r_h = (pred - y) / sigma
        if z is None:
            return r_h
        r_l = (load0 + v * T - z) / SIGMA_L
        return np.concatenate([r_h, r_l])

    lo = np.array([-40.0, 1e-5, 1e-4])
    hi = np.array([5.0, 1.0, 5.0])
    rows = []
    for factor in factors:
        gamma = float(GAMMA * factor)
        starts = [
            np.array([LOG_MU_B, VBAR, L0]),
            np.array([LOG_MU_B, VBAR * GAMMA / gamma, L0 * GAMMA / gamma]),
        ]
        best = None
        for x0 in starts:
            x0 = np.minimum(np.maximum(x0, lo + 1e-12), hi - 1e-12)
            fit = least_squares(
                residuals,
                x0,
                args=(gamma,),
                bounds=(lo, hi),
                method="trf",
                xtol=1e-14,
                ftol=1e-14,
                gtol=1e-14,
                max_nfev=80,
            )
            rss = float(np.dot(fit.fun, fit.fun))
            rec = {
                "gamma": gamma,
                "factor": float(factor),
                "rss": rss,
                "success": bool(fit.success),
                "log_mu_b": float(fit.x[0]),
                "v": float(fit.x[1]),
                "L0": float(fit.x[2]),
                "start": "truth" if np.allclose(x0, starts[0]) else "compensator",
            }
            if best is None or rss < best["rss"]:
                best = rec
        rows.append(best)
    baseline = min(row["rss"] for row in rows)
    for row in rows:
        row["delta_chi2"] = row["rss"] - baseline
    return rows


def profile_beta(y: np.ndarray, sigma: np.ndarray) -> list[dict]:
    weight = 1.0 / sigma**2
    rows = []
    beta_hat = float(wls_poly(T, y, sigma, 1)[0][1])
    for factor in BETA_FACTORS:
        beta = float(beta_hat * factor) if abs(beta_hat) > 0 else float(factor)
        # The plotted grid is relative to the generating slope, not the noisy MLE,
        # so a flat generating truth stays centred. Recompute against GAMMA*VBAR.
        beta = float(GAMMA * VBAR * factor)
        alpha = float(np.sum(weight * (y - beta * T)) / np.sum(weight))
        rss = float(np.sum(((y - alpha - beta * T) / sigma) ** 2))
        rows.append({"beta": beta, "factor": float(factor), "alpha": alpha, "rss": rss})
    baseline = min(row["rss"] for row in rows)
    for row in rows:
        row["delta_chi2"] = row["rss"] - baseline
    return rows


def _cross(x0: float, d0: float, x1: float, d1: float, level: float) -> float:
    if d1 == d0:
        return float(x1)
    weight = (level - d0) / (d1 - d0)
    return float(x0 + weight * (x1 - x0))


def profile_interval(rows: list[dict], key: str) -> dict:
    xs = np.array([row[key] for row in rows], dtype=float)
    ds = np.array([row["delta_chi2"] for row in rows], dtype=float)
    order = np.argsort(xs)
    xs, ds = xs[order], ds[order]
    imin = int(np.argmin(ds))
    lo = None
    hi = None
    for i in range(imin, 0, -1):
        if ds[i] <= CHI2_95_1 < ds[i - 1] or ds[i - 1] <= CHI2_95_1 < ds[i]:
            lo = _cross(float(xs[i]), float(ds[i]), float(xs[i - 1]), float(ds[i - 1]), CHI2_95_1)
            break
    for i in range(imin, xs.size - 1):
        if ds[i] <= CHI2_95_1 < ds[i + 1] or ds[i + 1] <= CHI2_95_1 < ds[i]:
            hi = _cross(float(xs[i]), float(ds[i]), float(xs[i + 1]), float(ds[i + 1]), CHI2_95_1)
            break
    closed = bool(ds[0] > CHI2_95_1 and ds[-1] > CHI2_95_1 and lo is not None and hi is not None)
    inside = ds <= CHI2_95_1
    return {
        "closed_on_grid": closed,
        "lo": lo,
        "hi": hi,
        "grid_inside_lo": None if not np.any(inside) else float(xs[inside].min()),
        "grid_inside_hi": None if not np.any(inside) else float(xs[inside].max()),
        "max_delta": float(ds.max()),
        "min_delta": float(ds.min()),
        "endpoint_delta": [float(ds[0]), float(ds[-1])],
    }


def profile_g(y: np.ndarray, sigma: np.ndarray, g_grid: np.ndarray) -> list[dict]:
    def residuals(x: np.ndarray, g: float) -> np.ndarray:
        log_mu_b, gamma, v, load0 = (float(v_) for v_ in x)
        pred = log_mu_b + gamma * mean_damage_scalar(T, v, load0, g)
        return (pred - y) / sigma

    lo = np.array([-40.0, 0.05, 1e-5, 1e-4])
    hi = np.array([5.0, 80.0, 1.0, 5.0])
    x_true = np.array([LOG_MU_B, GAMMA, VBAR, L0])
    rows = []
    for g in g_grid:
        fit = least_squares(
            residuals,
            x_true,
            args=(float(g),),
            bounds=(lo, hi),
            method="trf",
            xtol=1e-14,
            ftol=1e-14,
            gtol=1e-14,
            max_nfev=120,
        )
        rss = float(np.dot(fit.fun, fit.fun))
        rows.append(
            {
                "g": float(g),
                "rss": rss,
                "success": bool(fit.success),
                "log_mu_b": float(fit.x[0]),
                "gamma": float(fit.x[1]),
                "v": float(fit.x[2]),
                "L0": float(fit.x[3]),
            }
        )
    baseline = min(row["rss"] for row in rows)
    for row in rows:
        row["delta_chi2"] = row["rss"] - baseline
    return rows


def _dm_partials(t: np.ndarray, v: float, load0: float, g: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mean damage and partials with respect to v, L0 and g."""
    t = np.asarray(t, dtype=float)
    gt = g * t
    if abs(g) < 1e-8:
        # m = L0*(1 + gt + (gt)^2/2 + (gt)^3/6) + v*t*(1 + gt/2 + (gt)^2/6 + (gt)^3/24)
        m = mean_damage_scalar(t, v, load0, g)
        dm_dv = t * (1.0 + gt / 2.0 + gt**2 / 6.0 + gt**3 / 24.0 + gt**4 / 120.0)
        dm_dl = np.exp(gt)
        dm_dg = load0 * t * (1.0 + gt + gt**2 / 2.0 + gt**3 / 6.0) + v * (t**2) * (
            0.5 + gt / 3.0 + gt**2 / 8.0 + gt**3 / 30.0
        )
        return m, dm_dv, dm_dl, dm_dg
    expo = np.exp(gt)
    m = expo * load0 + v * np.expm1(gt) / g
    dm_dv = np.expm1(gt) / g
    dm_dl = expo
    dm_dg = t * expo * load0 + v * (t * expo * g - expo + 1.0) / g**2
    return m, dm_dv, dm_dl, dm_dg


def jac_g_model(g: float) -> tuple[np.ndarray, np.ndarray]:
    """Analytic Jacobian of log hazard for (log_mu_b, gamma, v, L0, g)."""
    m, dm_dv, dm_dl, dm_dg = _dm_partials(T, VBAR, L0, g)
    jac = np.column_stack(
        [
            np.ones_like(T),
            m,
            GAMMA * dm_dv,
            GAMMA * dm_dl,
            GAMMA * dm_dg,
        ]
    )
    mu = np.exp(LOG_MU_B + GAMMA * m)
    sigma = sigma_log_from_mu(mu, exposure(T))
    return jac, sigma


def fmt(x: float, digits: int = 6) -> float:
    return float(f"{x:.{digits}g}")


def main() -> None:
    expo = exposure(T)
    log_mu_true = log_mu_gain(T, LOG_MU_B, GAMMA, VBAR, L0)
    mu_true = mu_from_log(log_mu_true)
    sigma = sigma_log_from_mu(mu_true, expo)
    damage = damage_linear(T)
    load = damage.mean(axis=1)

    # --- regimes: same linear damage, two hazard maps; plus coupled damage ---
    mu_add_lo = float(mu_true[0])
    mu_add_hi = float(mu_true[-1])
    load_lo = float(load[0])
    load_hi = float(load[-1])
    mu_add = mu_add_lo + (mu_add_hi - mu_add_lo) * (load - load_lo) / (load_hi - load_lo)
    log_mu_add = np.log(mu_add)

    call_gain = gompertz_call(log_mu_true, sigma)
    call_add = gompertz_call(log_mu_add, sigma)
    lin_damage = node_linearity(damage)

    scan = []
    for g in G_SCAN:
        net = damage_network(T, float(g))
        mean_d = net.mean(axis=1)
        scalar = mean_damage_scalar(T, VBAR, L0, float(g))
        log_mu_g = LOG_MU_B + GAMMA * mean_d
        # Lack-of-fit uses the exposure weights of the *linear* generator, so the
        # noise budget is the demographic schedule under study, not a moving target.
        call = gompertz_call(log_mu_g, sigma)
        nodes = node_linearity(net)
        scan.append(
            {
                "g": float(g),
                "mean_matches_scalar": float(np.max(np.abs(mean_d - scalar))),
                "damage_min_r2": nodes["min_r2"],
                "damage_max_rel_resid": nodes["max_rel_resid"],
                "damage_near_linear": nodes["near_linear"],
                "gompertz_like": call["gompertz_like"],
                "rss_linear": call["rss_linear"],
                "lr_vs_quadratic": call["lr_vs_quadratic"],
                "quadratic_c": call["quadratic_c"],
                "beta_hat": call["beta"],
                "log10_mu_30": float(np.log10(np.exp(log_mu_g[0]))),
                "log10_mu_95": float(np.log10(np.exp(log_mu_g[-1]))),
            }
        )

    # --- graphs ---
    graphs = [
        graph_summary("uncoupled", np.zeros((N, N)), compensate=True),
        graph_summary("ring", ring_adjacency(0.45), compensate=True),
        graph_summary("star", star_adjacency(0.30), compensate=True),
        graph_summary("ring_gamma_fixed", ring_adjacency(0.45), compensate=False),
        graph_summary("star_gamma_fixed", star_adjacency(0.30), compensate=False),
    ]

    # --- Fisher, mechanistic scalar and node vector ---
    fishers = {}
    for label, with_load in (("P", False), ("PL", True)):
        jac, sig = jac_mechanistic_4(GAMMA, VBAR, L0, with_load)
        fisher = fisher_outer(jac, sig)
        fishers[label] = {
            "matrix_report": eig_report(fisher),
            "columns": column_spectrum(jac, sig),
            "cr_var": cr_diag(fisher),
        }
        if label == "P":
            n1 = np.array([-GAMMA, 0.0, 0.0, 1.0])
            n2 = np.array([0.0, 1.0, -VBAR / GAMMA, -L0 / GAMMA])
            fishers[label]["null_alignment"] = null_alignment(fisher, [n1, n2])
            fishers[label]["analytic_null"] = {
                "baseline_load": n1.tolist(),
                "gain_velocity": n2.tolist(),
            }
        if fishers[label]["cr_var"] is not None:
            se = np.sqrt(fishers[label]["cr_var"])
            theta = np.array([LOG_MU_B, GAMMA, VBAR, L0])
            fishers[label]["rel_se"] = (se / np.abs(theta)).tolist()
            fishers[label]["se"] = se.tolist()

    theta14 = np.concatenate([np.array([LOG_MU_B, GAMMA]), V, D0])
    for label, kind in (("P14", "hazard"), ("PL14", "mean"), ("PN14", "nodes")):
        jac, sig = jac_nodes(kind)
        fisher = fisher_outer(jac, sig)
        fishers[label] = {
            "matrix_report": eig_report(fisher),
            "columns": column_spectrum(jac, sig),
            "cr_var": cr_diag(fisher),
        }
        if fishers[label]["cr_var"] is not None:
            se = np.sqrt(fishers[label]["cr_var"])
            fishers[label]["rel_se"] = (se / np.abs(theta14)).tolist()

    # Observable (alpha, beta) Fisher. log mu = alpha + beta * t.
    j_obs = np.column_stack([np.ones_like(T), T])
    f_obs = fisher_outer(j_obs, sigma)
    obs_report = eig_report(f_obs)
    obs_columns = column_spectrum(j_obs, sigma)
    obs_var = cr_diag(f_obs)
    alpha_true = LOG_MU_B + GAMMA * L0
    beta_true = GAMMA * VBAR
    obs_se = np.sqrt(obs_var) if obs_var is not None else None

    # Coupled-model Fisher at the two pre-specified gains.
    coupled_fisher = {}
    for g in G_PROFILE_TRUTHS:
        jac, sig = jac_g_model(float(g))
        fisher = fisher_outer(jac, sig)
        report = eig_report(fisher)
        evals, evecs = np.linalg.eigh(0.5 * (fisher + fisher.T))
        null_g = evecs[4, evals < ABS_EIG_FLOOR]
        n1 = np.array([-GAMMA, 0.0, 0.0, 1.0, 0.0])
        n2 = np.array([0.0, 1.0, -VBAR / GAMMA, -L0 / GAMMA, 0.0])
        weighted = jac / sig[:, None]
        coupled_fisher[str(g)] = {
            **report,
            "columns": column_spectrum(jac, sig),
            "null_g_components": [float(v) for v in np.atleast_1d(null_g)],
            "tradeoff_weighted_residual": {
                "baseline_load": float(np.linalg.norm(weighted @ n1)),
                "gain_velocity": float(np.linalg.norm(weighted @ n2)),
            },
        }

    # --- profiles ---
    z_true = L0 + VBAR * T
    prof_g_P = profile_gamma(log_mu_true, sigma, None)
    prof_g_PL = profile_gamma(log_mu_true, sigma, z_true)
    prof_g_PL_fine = profile_gamma(log_mu_true, sigma, z_true, factors=GAMMA_FINE_FACTORS)
    prof_b = profile_beta(log_mu_true, sigma)

    rng = np.random.default_rng(SEED)
    y_noisy = log_mu_true + rng.normal(0.0, sigma)
    z_noisy = z_true + rng.normal(0.0, SIGMA_L, size=T.size)
    prof_g_P_noisy = profile_gamma(y_noisy, sigma, None)
    prof_b_noisy = profile_beta(y_noisy, sigma)
    prof_g_PL_noisy = profile_gamma(y_noisy, sigma, z_noisy)
    prof_g_PL_noisy_fine = profile_gamma(y_noisy, sigma, z_noisy, factors=GAMMA_FINE_FACTORS)

    g_grid = np.linspace(-0.02, 0.08, 41)
    g_profiles = {}
    g_profiles_fine = {}
    for g_true in G_PROFILE_TRUTHS:
        y_g = LOG_MU_B + GAMMA * mean_damage_scalar(T, VBAR, L0, float(g_true))
        # Noise budget stays the linear-generator sigma, as in the scan.
        g_profiles[str(g_true)] = profile_g(y_g, sigma, g_grid)
        g_profiles_fine[str(g_true)] = profile_g(y_g, sigma, G_FINE[float(g_true)])

    # Gamma profile along the ridge should keep the product gamma*v at beta.
    products = [row["gamma"] * row["v"] for row in prof_g_P]

    def pack_profile(rows: list[dict], key: str) -> dict:
        return {
            "grid": [{k: row[k] for k in (key, "factor", "delta_chi2", "rss") if k in row} | (
                {"g": row["g"]} if "g" in row and key != "g" else {}
            ) for row in rows] if False else [
                {k: row[k] for k in row if k in (key, "factor", "delta_chi2", "rss", "g", "gamma", "beta", "v", "L0")}
                for row in rows
            ],
            "interval": profile_interval(rows, key),
            "max_delta": float(max(row["delta_chi2"] for row in rows)),
        }

    # Truth-start RSS on the noiseless demographic gamma profile (second start is recorded
    # only as the winner). Recompute the truth-start RSS explicitly for the assertion.
    truth_start_rss = []
    lo = np.array([-40.0, 1e-5, 1e-4])
    hi = np.array([5.0, 1.0, 5.0])

    def resid_truth(x: np.ndarray, gamma: float) -> np.ndarray:
        return (log_mu_gain(T, float(x[0]), gamma, float(x[1]), float(x[2])) - log_mu_true) / sigma

    for factor in (0.5, 2.0):
        fit = least_squares(
            resid_truth,
            np.array([LOG_MU_B, VBAR, L0]),
            args=(float(GAMMA * factor),),
            bounds=(lo, hi),
            method="trf",
            xtol=1e-14,
            ftol=1e-14,
            gtol=1e-14,
            max_nfev=80,
        )
        truth_start_rss.append(float(np.dot(fit.fun, fit.fun)))

    results = {
        "seed": SEED,
        "n_nodes": N,
        "ages": [float(A0), float(A1), int(AGES.size)],
        "truth": {
            "v": V.tolist(),
            "D0": D0.tolist(),
            "v_bar": VBAR,
            "L0": L0,
            "gamma": GAMMA,
            "mu_b": MU_B,
            "log_mu_b": LOG_MU_B,
            "alpha": alpha_true,
            "beta": beta_true,
            "doubling_years": float(np.log(2.0) / beta_true),
            "mu_age_30": float(mu_true[0]),
            "mu_age_60": float(mu_true[int(60 - A0)]),
            "mu_age_95": float(mu_true[-1]),
            "sigma_log_age_30": float(sigma[0]),
            "sigma_log_age_60": float(sigma[int(60 - A0)]),
            "sigma_log_age_95": float(sigma[-1]),
            "exposure_age_30": float(expo[0]),
            "exposure_age_95": float(expo[-1]),
            "fold_mu_95_over_30": float(mu_true[-1] / mu_true[0]),
        },
        "cuts": {
        "abs_eig_floor": ABS_EIG_FLOOR,
        "column_sv_tol": COLUMN_SV_TOL,
            "chi2_95_1": CHI2_95_1,
            "flat_chi2": FLAT_CHI2,
            "near_linear_r2": NEAR_LINEAR_R2,
            "near_linear_max_rel": NEAR_LINEAR_MAX_REL,
            "sigma_L": SIGMA_L,
            "sigma_D": SIGMA_D,
        },
        "regimes": {
            "linear_damage": lin_damage,
            "load_x_gain": call_gain,
            "additive_endpoint_matched": call_add,
        },
        "g_scan": scan,
        "graphs": [
            {k: v for k, v in g.items() if k not in ("log_mu", "felt")}
            for g in graphs
        ],
        "fisher": fishers,
        "observable_gompertz": {
            "report": obs_report,
            "columns": obs_columns,
            "se_alpha": None if obs_se is None else float(obs_se[0]),
            "se_beta": None if obs_se is None else float(obs_se[1]),
            "rel_se_beta": None if obs_se is None else float(obs_se[1] / abs(beta_true)),
            "profile_ci_beta": profile_interval(prof_b, "beta"),
        },
        "coupled_fisher": coupled_fisher,
        "profiles": {
            "gamma_P": pack_profile(prof_g_P, "gamma"),
            "gamma_PL": pack_profile(prof_g_PL, "gamma"),
            "gamma_PL_fine": pack_profile(prof_g_PL_fine, "gamma"),
            "beta_P": pack_profile(prof_b, "beta"),
            "gamma_P_noisy": pack_profile(prof_g_P_noisy, "gamma"),
            "gamma_PL_noisy": pack_profile(prof_g_PL_noisy, "gamma"),
            "gamma_PL_noisy_fine": pack_profile(prof_g_PL_noisy_fine, "gamma"),
            "beta_P_noisy": pack_profile(prof_b_noisy, "beta"),
            "gamma_times_v_on_P": products,
            "truth_start_rss_half_and_double": truth_start_rss,
            "g_ode": {key: pack_profile(rows, "g") for key, rows in g_profiles.items()},
            "g_ode_fine": {key: pack_profile(rows, "g") for key, rows in g_profiles_fine.items()},
        },
    }

    # ---------- figures ----------
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 140,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
        }
    )
    colours = {
        "gain": "#1f4e79",
        "add": "#b85c38",
        "mild": "#2f6f4e",
        "strong": "#6b4c9a",
    }

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5))
    axes[0].plot(AGES, load, color=colours["gain"], lw=1.8, label="Mean load, uncoupled")
    net_mild = damage_network(T, 0.005).mean(axis=1)
    net_strong = damage_network(T, 0.04).mean(axis=1)
    axes[0].plot(AGES, net_mild, color=colours["mild"], lw=1.4, label="Mean load, g = 0.005")
    axes[0].plot(AGES, net_strong, color=colours["strong"], lw=1.4, label="Mean load, g = 0.04")
    axes[0].set_xlabel("Age (model years)")
    axes[0].set_ylabel("Mean damage load")
    axes[0].legend(frameon=False)
    axes[1].plot(AGES, np.log10(mu_true), color=colours["gain"], lw=1.8, label="Load × gain")
    axes[1].plot(AGES, np.log10(mu_add), color=colours["add"], lw=1.4, label="Additive, ends matched")
    axes[1].plot(AGES, np.log10(np.exp(LOG_MU_B + GAMMA * net_mild)), color=colours["mild"], lw=1.4, label="g = 0.005")
    axes[1].plot(AGES, np.log10(np.exp(LOG_MU_B + GAMMA * net_strong)), color=colours["strong"], lw=1.4, label="g = 0.04")
    axes[1].set_xlabel("Age (model years)")
    axes[1].set_ylabel("log10 hazard")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR / "damage_and_hazard.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5))
    gs = np.array([row["g"] for row in scan])
    axes[0].plot(gs, [row["damage_max_rel_resid"] for row in scan], "o-", color=colours["gain"], lw=1.4)
    axes[0].axhline(NEAR_LINEAR_MAX_REL, color="#888", ls="--", lw=0.8)
    axes[0].set_xlabel("Damage-rate gain g")
    axes[0].set_ylabel("Max relative linear residual")
    axes[1].plot(gs, [row["rss_linear"] for row in scan], "o-", color=colours["add"], lw=1.4)
    axes[1].axhline(call_gain["chi2_crit_df"], color="#888", ls="--", lw=0.8)
    axes[1].set_yscale("symlog", linthresh=1.0)
    axes[1].set_xlabel("Damage-rate gain g")
    axes[1].set_ylabel("Gompertz weighted RSS")
    fig.tight_layout()
    fig.savefig(FIGDIR / "coupling_scan.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    series = [
        ("Hazard only, 4 params", fishers["P"]["columns"]["singular_values"]),
        ("Hazard + mean load, 4 params", fishers["PL"]["columns"]["singular_values"]),
        ("Hazard only, 14 params", fishers["P14"]["columns"]["singular_values"]),
        ("Hazard + node damage, 14 params", fishers["PN14"]["columns"]["singular_values"]),
    ]
    for label, sv in series:
        sv = np.array(sv, dtype=float)
        if sv.size == 0:
            continue
        ax.plot(np.arange(1, sv.size + 1), sv / sv[0], "o-", ms=3.5, lw=1.1, label=label)
    ax.axhline(COLUMN_SV_TOL, color="#888", ls="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Singular-value index")
    ax.set_ylabel("Column-balanced singular value / largest")
    ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fisher_spectra.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.2))
    axes[0, 0].plot(
        [row["factor"] for row in prof_g_P],
        [row["delta_chi2"] for row in prof_g_P],
        "o-",
        color=colours["add"],
        lw=1.4,
    )
    axes[0, 0].axhline(CHI2_95_1, color="#888", ls="--", lw=0.8)
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_ylim(-0.05, 1.0)
    axes[0, 0].set_xlabel("γ / γ true, hazard only")
    axes[0, 0].set_ylabel("Profile Δχ²")
    axes[0, 1].plot(
        [row["factor"] for row in prof_g_PL_fine],
        [row["delta_chi2"] for row in prof_g_PL_fine],
        "o-",
        color=colours["gain"],
        lw=1.4,
        ms=3.5,
    )
    axes[0, 1].axhline(CHI2_95_1, color="#888", ls="--", lw=0.8)
    axes[0, 1].set_xlabel("γ / γ true, hazard + mean load")
    axes[0, 1].set_ylabel("Profile Δχ²")
    axes[1, 0].plot(
        [row["factor"] for row in prof_b],
        [row["delta_chi2"] for row in prof_b],
        "o-",
        color=colours["gain"],
        lw=1.4,
        ms=3.5,
    )
    axes[1, 0].axhline(CHI2_95_1, color="#888", ls="--", lw=0.8)
    axes[1, 0].set_xlabel("β / β true, hazard only")
    axes[1, 0].set_ylabel("Profile Δχ²")
    axes[1, 1].plot(
        [row["g"] for row in g_profiles_fine["0.0"]],
        [row["delta_chi2"] for row in g_profiles_fine["0.0"]],
        "o-",
        color=colours["strong"],
        lw=1.4,
        ms=3.5,
    )
    axes[1, 1].axhline(CHI2_95_1, color="#888", ls="--", lw=0.8)
    axes[1, 1].set_xlabel("Damage-rate gain g, truth 0")
    axes[1, 1].set_ylabel("Profile Δχ²")
    fig.tight_layout()
    fig.savefig(FIGDIR / "profiles.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5))
    for ginfo, colour, label in (
        (graphs[0], "#444", "Uncoupled"),
        (graphs[1], colours["gain"], "Ring, γ retuned"),
        (graphs[2], colours["add"], "Star, γ retuned"),
    ):
        axes[0].plot(AGES, ginfo["log_mu"] / np.log(10.0), color=colour, lw=1.5, label=label)
    axes[0].set_xlabel("Age (model years)")
    axes[0].set_ylabel("log10 hazard")
    axes[0].legend(frameon=False)
    # Node 0 felt load, compensated graphs still have different felt trajectories.
    for ginfo, colour, label in (
        (graphs[0], "#444", "Uncoupled"),
        (graphs[1], colours["gain"], "Ring"),
        (graphs[2], colours["add"], "Star"),
    ):
        axes[1].plot(AGES, ginfo["felt"][:, 0], color=colour, lw=1.5, label=label)
    axes[1].set_xlabel("Age (model years)")
    axes[1].set_ylabel("Felt load at node 1")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR / "graphs_same_hazard.png")
    plt.close(fig)

    # Strip bulky per-grid vectors? Keep them; the thesis cites aggregates.
    # Drop per-row start bookkeeping already not stored.
    out = ROOT / "results.json"
    # Profiles currently include full grids. Fine.
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")

    assertions = []

    def check(name: str, ok: bool) -> None:
        assertions.append({"name": name, "ok": bool(ok)})
        if not ok:
            raise AssertionError(name)

    check("damage near linear", lin_damage["near_linear"])
    check("gain is gompertz-like", call_gain["gompertz_like"])
    check("additive is not", not call_add["gompertz_like"])
    check("P rank 2 of 4", fishers["P"]["matrix_report"]["numerical_rank"] == 2)
    check("P column rank 2", fishers["P"]["columns"]["column_rank"] == 2)
    check("PL rank 4 of 4", fishers["PL"]["matrix_report"]["numerical_rank"] == 4)
    check("PL column rank 4", fishers["PL"]["columns"]["column_rank"] == 4)
    check("P14 rank 2", fishers["P14"]["matrix_report"]["numerical_rank"] == 2)
    check("P14 column rank 2", fishers["P14"]["columns"]["column_rank"] == 2)
    check("PL14 rank 4", fishers["PL14"]["matrix_report"]["numerical_rank"] == 4)
    check("PN14 rank 14", fishers["PN14"]["matrix_report"]["numerical_rank"] == 14)
    check("PN14 column rank 14", fishers["PN14"]["columns"]["column_rank"] == 14)
    check("coupled g=0 rank 3", coupled_fisher["0.0"]["numerical_rank"] == 3)
    check("coupled g=0.02 rank 3", coupled_fisher["0.02"]["numerical_rank"] == 3)
    check("null alignment", min(fishers["P"]["null_alignment"]) > 0.99)
    check("gamma profile flat", results["profiles"]["gamma_P"]["max_delta"] < FLAT_CHI2)
    check("truth start finds ridge", max(truth_start_rss) < 1e-6)
    check("beta profile closed", results["profiles"]["beta_P"]["interval"]["closed_on_grid"])
    check("gamma with load closed", results["profiles"]["gamma_PL"]["interval"]["closed_on_grid"])
    check("gamma fine ci", results["profiles"]["gamma_PL_fine"]["interval"]["closed_on_grid"])
    check("g fine ci", results["profiles"]["g_ode_fine"]["0.0"]["interval"]["closed_on_grid"])
    check("ring compensated", graphs[1]["max_abs_log_mu_gap"] < 1e-8)
    check("star compensated", graphs[2]["max_abs_log_mu_gap"] < 1e-8)
    check("ring fixed gamma visible", graphs[3]["max_abs_log_mu_gap"] > 0.5)
    check("scalar mean matches network", max(row["mean_matches_scalar"] for row in scan) < 1e-8)
    check("observable rank 2", obs_report["numerical_rank"] == 2)

    results["assertions"] = assertions
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print("beta", beta_true, "doubling", np.log(2) / beta_true)
    print("mu30/60/95", mu_true[0], mu_true[int(60 - A0)], mu_true[-1], "fold", mu_true[-1] / mu_true[0])
    print("sigma30/95", sigma[0], sigma[-1])
    print("gain call", {k: call_gain[k] for k in ("rss_linear", "beta", "gompertz_like", "lr_vs_quadratic")})
    print("add call", {k: call_add[k] for k in ("rss_linear", "beta", "gompertz_like", "quadratic_c", "lr_vs_quadratic")})
    print("damage", lin_damage["min_r2"], lin_damage["max_rel_resid"])
    print("SCAN")
    for row in scan:
        print(
            f"  g={row['g']:.3f} near={row['damage_near_linear']} gomp={row['gompertz_like']} "
            f"rel={row['damage_max_rel_resid']:.4g} rss={row['rss_linear']:.4g} "
            f"lr={row['lr_vs_quadratic']:.4g} c={row['quadratic_c']:.4g}"
        )
    print("FISHER")
    for key in ("P", "PL", "P14", "PL14", "PN14"):
        rep = fishers[key]["matrix_report"]
        cols = fishers[key]["columns"]["column_rank"]
        print(f"  {key} num={rep['numerical_rank']} col={cols} of {rep['n_param']} cond={rep['condition_nonzero']}")
    print("  null align", fishers["P"]["null_alignment"])
    if fishers["PL"]["cr_var"] is not None:
        print("  PL rel se", fishers["PL"]["rel_se"])
    print("  obs se beta", None if obs_se is None else float(obs_se[1]), "rel", None if obs_se is None else float(obs_se[1] / beta_true))
    print("  beta interval", results["observable_gompertz"]["profile_ci_beta"])
    print("GAMMA max dchi P/PL", results["profiles"]["gamma_P"]["max_delta"], results["profiles"]["gamma_PL"]["max_delta"])
    print("  PL interval", results["profiles"]["gamma_PL"]["interval"])
    print("  PL fine", results["profiles"]["gamma_PL_fine"]["interval"])
    print("  g fine", {k: v["interval"] for k, v in results["profiles"]["g_ode_fine"].items()})
    print("  noisy gamma max", results["profiles"]["gamma_P_noisy"]["max_delta"], "noisy beta", results["profiles"]["beta_P_noisy"]["interval"])
    print("  noisy PL", results["profiles"]["gamma_PL_noisy"]["interval"])
    print("G ODE FISHER")
    for key, rep in coupled_fisher.items():
        print(
            f"  g={key} num={rep['numerical_rank']} col={rep['columns']['column_rank']} "
            f"of {rep['n_param']} null_g={['%.3g' % v for v in rep['null_g_components']]} "
            f"resid={rep['tradeoff_weighted_residual']}"
        )
    for key, packed in results["profiles"]["g_ode"].items():
        print(f"  profile g truth {key}", packed["interval"])
    print("GRAPHS")
    for ginfo in graphs:
        print(f"  {ginfo['name']} gamma={ginfo['gamma']:.6g} gap={ginfo['max_abs_log_mu_gap']:.3g} slope={ginfo['load_slope']:.6g}")
    print("products gamma*v", min(products), max(products), "target", beta_true)
    print("truth start rss", truth_start_rss)
    print("assertions", len(assertions), "ok")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
