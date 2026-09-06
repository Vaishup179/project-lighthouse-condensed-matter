from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import optimize, signal


DATA = Path("/app/data")
OUT = Path("/app/output")
ALPHA_LK = 14.694


def fit_line(path: Path, sensor: str, reference: str) -> tuple[float, float]:
    frame = pd.read_csv(path)
    gain, offset = np.polyfit(frame[sensor].to_numpy(float), frame[reference].to_numpy(float), 1)
    return float(gain), float(offset)


@pytest.fixture(scope="session")
def calibrations() -> dict[str, float]:
    fg, fo = fit_line(DATA / "field_calibration.csv", "sensor_field_t", "reference_field_t")
    ag, ao = fit_line(DATA / "angle_calibration.csv", "encoder_angle_deg", "reference_angle_deg")
    tg, to = fit_line(DATA / "temperature_calibration.csv", "sensor_temperature_k", "reference_temperature_k")
    return {"fg": fg, "fo": fo, "ag": ag, "ao": ao, "tg": tg, "to": to}


def verifier_trace(group: pd.DataFrame, cal: dict[str, float]) -> tuple[np.ndarray, np.ndarray, float, float]:
    B = cal["fg"] * group["field_sensor_t"].to_numpy(float) + cal["fo"]
    angle = cal["ag"] * group["encoder_deg"].to_numpy(float) + cal["ao"]
    temp = cal["tg"] * group["temperature_sensor_k"].to_numpy(float) + cal["to"]
    y = group["rho_xx_uohm_cm"].to_numpy(float)

    good = np.isfinite(B) & np.isfinite(y)
    B, y = B[good], y[good]
    order = np.argsort(B)
    B, y = B[order], y[order]

    # Collapse exact duplicate field readings without using any oracle code.
    rounded = np.round(B, 10)
    unique, inverse = np.unique(rounded, return_inverse=True)
    B2 = np.zeros(unique.size)
    y2 = np.zeros(unique.size)
    for i in range(unique.size):
        use = inverse == i
        B2[i] = np.median(B[use])
        y2[i] = np.median(y[use])
    B, y = B2, y2

    keep = np.ones(B.size, dtype=bool)
    for _ in range(4):
        coeff = np.polyfit(B[keep], y[keep], 4)
        residual = y - np.polyval(coeff, B)
        scale = 1.4826 * np.median(np.abs(residual[keep] - np.median(residual[keep])))
        if scale <= 0:
            break
        new_keep = np.abs(residual - np.median(residual[keep])) < 6.5 * scale
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    residual = y - np.polyval(coeff, B)

    # Hampel-like local replacement for isolated acquisition spikes.
    med = signal.medfilt(residual, kernel_size=7)
    delta = residual - med
    local_scale = 1.4826 * np.median(np.abs(delta - np.median(delta)))
    if local_scale > 0:
        bad = np.abs(delta - np.median(delta)) > 8.0 * local_scale
        residual[bad] = med[bad]

    x = 1.0 / B
    order = np.argsort(x)
    x, residual = x[order], residual[order]
    xu = np.linspace(x.min(), x.max(), 3072)
    yu = np.interp(xu, x, residual)
    yu = signal.detrend(yu)
    return xu, yu, float(np.nanmedian(angle)), float(np.nanmedian(temp))


def verifier_frequency(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    window = np.hanning(x.size)
    yw = y * window
    freqs = np.fft.rfftfreq(x.size, d=x[1] - x[0])
    amp = 2.0 * np.abs(np.fft.rfft(yw)) / window.sum()
    band = (freqs >= 80.0) & (freqs <= 900.0)
    idx = np.flatnonzero(band)[np.argmax(amp[band])]
    coarse = freqs[idx]
    df = freqs[1] - freqs[0]

    def objective(freq: float) -> float:
        z = np.sum(yw * np.exp(-2j * np.pi * freq * x))
        return float(-np.abs(z))

    refined = optimize.minimize_scalar(
        objective,
        bounds=(max(80.0, coarse - df), min(900.0, coarse + df)),
        method="bounded",
    )
    frequency = float(refined.x)
    coeff = np.sum(yw * np.exp(-2j * np.pi * frequency * x))
    amplitude = float(2.0 * np.abs(coeff) / window.sum())
    return frequency, amplitude


@pytest.fixture(scope="session")
def raw_reconstruction(calibrations: dict[str, float]):
    sweeps = pd.read_csv(DATA / "oscillation_sweeps.csv")
    trace = {}
    rows = []
    for sweep_id, group in sweeps.groupby("sweep_id", sort=True):
        x, y, angle, temp = verifier_trace(group, calibrations)
        freq, amp = verifier_frequency(x, y)
        trace[sweep_id] = (x, y)
        rows.append(
            {
                "sweep_id": sweep_id,
                "series": str(group["series"].iloc[0]),
                "angle": angle,
                "temp": temp,
                "freq": freq,
                "amp": amp,
            }
        )
    return pd.DataFrame(rows), trace


@pytest.fixture(scope="session")
def expected_angle(raw_reconstruction):
    rows, _ = raw_reconstruction
    angle = rows.loc[rows["series"] == "angle"].copy()
    angle["cluster"] = np.rint(angle["angle"]).astype(int)
    grouped = (
        angle.groupby("cluster", as_index=False)
        .agg(angle_deg=("angle", "mean"), frequency_t=("freq", "mean"))
        .sort_values("angle_deg")
        .reset_index(drop=True)
    )
    return grouped


def model_fit(expected_angle: pd.DataFrame) -> tuple[str, float]:
    theta = np.deg2rad(expected_angle["angle_deg"].to_numpy(float))
    observed = expected_angle["frequency_t"].to_numpy(float)

    q = 1.0 / np.abs(np.cos(theta))
    q_f0 = float(np.dot(q, observed) / np.dot(q, q))
    q_rss = float(np.sum((observed - q_f0 * q) ** 2))

    gamma = 3.2
    e = 1.0 / np.sqrt(np.cos(theta) ** 2 + np.sin(theta) ** 2 / gamma**2)
    e_f0 = float(np.dot(e, observed) / np.dot(e, e))
    e_rss = float(np.sum((observed - e_f0 * e) ** 2))

    if q_rss <= e_rss:
        return "quasi-2d", q_f0
    return "ellipsoidal-3d", e_f0


def verifier_segment_points(raw_reconstruction, expected_f0: float) -> pd.DataFrame:
    rows, traces = raw_reconstruction
    zero = rows[
        (rows["series"] == "temperature")
        | ((rows["series"] == "angle") & (np.abs(rows["angle"]) < 1.0))
    ]
    points = []
    for _, row in zero.iterrows():
        x, y = traces[str(row["sweep_id"])]
        left = int(0.07 * x.size)
        right = int(0.93 * x.size)
        edges = np.linspace(left, right, 7, dtype=int)  # six independent windows
        for seg, (start, stop) in enumerate(zip(edges[:-1], edges[1:])):
            xs, ys = x[start:stop], y[start:stop]
            w = np.hanning(xs.size)
            z = np.sum(ys * w * np.exp(-2j * np.pi * expected_f0 * xs))
            points.append(
                {
                    "temp": float(row["temp"]),
                    "temp_group": int(round(float(row["temp"]) * 10.0)),
                    "seg": seg,
                    "B": float(1.0 / np.mean(xs)),
                    "A": float(2.0 * np.abs(z) / w.sum()),
                }
            )
    frame = pd.DataFrame(points)
    return (
        frame.groupby(["temp_group", "seg"], as_index=False)
        .agg(temp=("temp", "mean"), B=("B", "mean"), A=("A", "mean"))
    )


def fit_mass_dingle(points: pd.DataFrame) -> tuple[float, float]:
    B = points["B"].to_numpy(float)
    T = points["temp"].to_numpy(float)
    A = points["A"].to_numpy(float)

    def residual(params: np.ndarray) -> np.ndarray:
        scale, mass, td = params
        x = ALPHA_LK * mass * T / B
        thermal = x / np.sinh(x)
        pred = scale * thermal * np.exp(-ALPHA_LK * mass * td / B)
        return np.log(np.maximum(A, 1e-14)) - np.log(np.maximum(pred, 1e-14))

    fit = optimize.least_squares(
        residual,
        x0=np.array([0.1, 0.30, 3.0]),
        bounds=(
            np.array([0.01, 0.05, 0.1]),
            np.array([1.0, 2.0, 20.0]),
        ),
    )
    return float(fit.x[1]), float(fit.x[2])


def load_summary() -> dict:
    return json.loads((OUT / "pocket_summary.json").read_text(encoding="utf-8"))


def test_required_outputs_and_schemas():
    required = [
        OUT / "angle_frequencies.csv",
        OUT / "holdout_predictions.csv",
        OUT / "pocket_summary.json",
        OUT / "analysis_report.md",
    ]
    for path in required:
        assert path.is_file(), f"Missing required output: {path}"

    angle = pd.read_csv(OUT / "angle_frequencies.csv")
    assert list(angle.columns) == [
        "angle_deg", "frequency_t", "frequency_spread_t", "peak_snr"
    ]
    assert len(angle) == 7

    holdout = pd.read_csv(OUT / "holdout_predictions.csv")
    assert list(holdout.columns) == ["angle_deg", "predicted_frequency_t"]
    assert len(holdout) == 3

    summary = load_summary()
    required_keys = {
        "field_gain", "field_offset_t", "angle_gain", "angle_offset_deg",
        "temperature_gain", "temperature_offset_k", "model", "F0_t",
        "effective_mass_me", "dingle_temperature_k", "extremal_area_m_inv2",
        "kF_m_inv", "quantum_lifetime_s", "fermi_velocity_m_s",
        "mean_free_path_m",
    }
    assert set(summary) == required_keys


def test_instrument_calibrations_are_recovered(calibrations):
    s = load_summary()
    assert np.isclose(s["field_gain"], calibrations["fg"], rtol=0, atol=2e-4)
    assert np.isclose(s["field_offset_t"], calibrations["fo"], rtol=0, atol=2e-3)
    assert np.isclose(s["angle_gain"], calibrations["ag"], rtol=0, atol=2e-4)
    assert np.isclose(s["angle_offset_deg"], calibrations["ao"], rtol=0, atol=0.04)
    assert np.isclose(s["temperature_gain"], calibrations["tg"], rtol=0, atol=3e-4)
    assert np.isclose(s["temperature_offset_k"], calibrations["to"], rtol=0, atol=0.004)


def test_angle_frequencies_match_raw_sweeps(expected_angle):
    submitted = pd.read_csv(OUT / "angle_frequencies.csv").sort_values("angle_deg").reset_index(drop=True)
    assert np.allclose(submitted["angle_deg"], expected_angle["angle_deg"], atol=0.35, rtol=0)
    assert np.allclose(submitted["frequency_t"], expected_angle["frequency_t"], atol=4.0, rtol=0)
    assert np.all(np.isfinite(submitted["frequency_spread_t"]))
    assert np.all(submitted["frequency_spread_t"] >= 0)
    assert np.all(np.isfinite(submitted["peak_snr"]))
    assert np.all(submitted["peak_snr"] > 0)


def test_model_and_zero_angle_frequency_are_supported(expected_angle):
    expected_model, expected_f0 = model_fit(expected_angle)
    s = load_summary()
    assert s["model"] == expected_model
    assert np.isclose(s["F0_t"], expected_f0, atol=3.5, rtol=0)


def test_mass_and_dingle_temperature_from_raw_damping(raw_reconstruction, expected_angle):
    _, expected_f0 = model_fit(expected_angle)
    points = verifier_segment_points(raw_reconstruction, expected_f0)
    expected_mass, expected_td = fit_mass_dingle(points)

    s = load_summary()
    assert np.isclose(s["effective_mass_me"], expected_mass, atol=0.055, rtol=0)
    assert np.isclose(s["dingle_temperature_k"], expected_td, atol=0.9, rtol=0)


def test_derived_quantities_are_physically_consistent():
    s = load_summary()
    c = json.loads((DATA / "constants.json").read_text(encoding="utf-8"))
    e = float(c["elementary_charge_c"])
    hbar = float(c["hbar_j_s"])
    kB = float(c["boltzmann_j_k"])
    me = float(c["electron_mass_kg"])

    area = 2.0 * np.pi * e * float(s["F0_t"]) / hbar
    kf = np.sqrt(area / np.pi)
    tau = hbar / (2.0 * np.pi * kB * float(s["dingle_temperature_k"]))
    vf = hbar * kf / (float(s["effective_mass_me"]) * me)
    mfp = vf * tau

    assert np.isclose(s["extremal_area_m_inv2"], area, rtol=0.015)
    assert np.isclose(s["kF_m_inv"], kf, rtol=0.015)
    assert np.isclose(s["quantum_lifetime_s"], tau, rtol=0.02)
    assert np.isclose(s["fermi_velocity_m_s"], vf, rtol=0.02)
    assert np.isclose(s["mean_free_path_m"], mfp, rtol=0.03)


def test_holdout_predictions_follow_selected_model():
    s = load_summary()
    submitted = pd.read_csv(OUT / "holdout_predictions.csv").sort_values("angle_deg").reset_index(drop=True)
    requested = pd.read_csv(DATA / "holdout_angles.csv").sort_values("angle_deg").reset_index(drop=True)

    assert np.allclose(submitted["angle_deg"], requested["angle_deg"], atol=1e-8, rtol=0)
    theta = np.deg2rad(requested["angle_deg"].to_numpy(float))
    F0 = float(s["F0_t"])

    if s["model"] == "quasi-2d":
        expected = F0 / np.abs(np.cos(theta))
    elif s["model"] == "ellipsoidal-3d":
        gamma = 3.2
        expected = F0 / np.sqrt(np.cos(theta)**2 + np.sin(theta)**2 / gamma**2)
    else:
        raise AssertionError(f"Unknown model: {s['model']}")

    assert np.allclose(submitted["predicted_frequency_t"], expected, rtol=0.01, atol=0.5)


def test_report_contains_a_scientific_summary():
    text = (OUT / "analysis_report.md").read_text(encoding="utf-8").strip()
    assert len(text) >= 350
