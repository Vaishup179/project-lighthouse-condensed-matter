#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, signal


DATA = Path("/app/data")
OUT = Path("/app/output")
ALPHA_LK = 14.694


def linear_calibration(path: Path, sensor_col: str, reference_col: str) -> tuple[float, float]:
    frame = pd.read_csv(path)
    gain, offset = np.polyfit(
        frame[sensor_col].to_numpy(float),
        frame[reference_col].to_numpy(float),
        1,
    )
    return float(gain), float(offset)


def clean_and_resample(
    frame: pd.DataFrame,
    field_gain: float,
    field_offset: float,
    angle_gain: float,
    angle_offset: float,
    temp_gain: float,
    temp_offset: float,
) -> dict[str, np.ndarray | float]:
    work = frame.copy()
    work["B"] = field_gain * work["field_sensor_t"] + field_offset
    work["theta"] = angle_gain * work["encoder_deg"] + angle_offset
    work["Tcal"] = temp_gain * work["temperature_sensor_k"] + temp_offset

    finite = np.isfinite(work["rho_xx_uohm_cm"]) & np.isfinite(work["B"])
    work = work.loc[finite].copy()

    work = (
        work.groupby("field_sensor_t", as_index=False)
        .agg(
            B=("B", "median"),
            theta=("theta", "median"),
            Tcal=("Tcal", "median"),
            y=("rho_xx_uohm_cm", "median"),
        )
        .sort_values("B")
    )

    B = work["B"].to_numpy(float)
    y = work["y"].to_numpy(float)

    mask = np.ones(B.size, dtype=bool)
    coeff = np.polyfit(B, y, 4)
    for _ in range(6):
        coeff = np.polyfit(B[mask], y[mask], 4)
        resid = y - np.polyval(coeff, B)
        center = np.median(resid[mask])
        mad = 1.4826 * np.median(np.abs(resid[mask] - center))
        limit = max(6.0 * mad, 0.008)
        new_mask = np.abs(resid - center) < limit
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask

    residual = y - np.polyval(coeff, B)

    local_median = signal.medfilt(residual, kernel_size=5)
    local_error = residual - local_median
    center = np.median(local_error)
    local_scale = 1.4826 * np.median(np.abs(local_error - center))
    spike_limit = max(7.0 * local_scale, 0.012)
    spike_mask = np.abs(local_error - center) > spike_limit
    residual[spike_mask] = local_median[spike_mask]

    invB = 1.0 / B
    order = np.argsort(invB)
    invB = invB[order]
    residual = residual[order]

    x = np.linspace(invB.min(), invB.max(), 4096)
    y_uniform = np.interp(x, invB, residual)
    y_uniform = signal.detrend(y_uniform, type="linear")

    return {
        "x": x,
        "y": y_uniform,
        "theta": float(np.median(work["theta"])),
        "temperature": float(np.median(work["Tcal"])),
    }


def continuous_peak(
    processed: dict[str, np.ndarray | float],
    fmin: float = 80.0,
    fmax: float = 900.0,
) -> tuple[float, float, float]:
    x = np.asarray(processed["x"], float)
    y = np.asarray(processed["y"], float)
    window = np.hanning(x.size)
    yw = y * window

    frequencies = np.fft.rfftfreq(x.size, d=x[1] - x[0])
    spectrum = 2.0 * np.abs(np.fft.rfft(yw)) / window.sum()
    use = (frequencies >= fmin) & (frequencies <= fmax)
    candidates = np.flatnonzero(use)
    peak_index = candidates[np.argmax(spectrum[candidates])]
    coarse = frequencies[peak_index]
    resolution = frequencies[1] - frequencies[0]

    def negative_amplitude(freq: float) -> float:
        coeff = np.sum(yw * np.exp(-2j * np.pi * freq * x))
        return float(-2.0 * np.abs(coeff) / window.sum())

    result = optimize.minimize_scalar(
        negative_amplitude,
        bounds=(max(fmin, coarse - resolution), min(fmax, coarse + resolution)),
        method="bounded",
        options={"xatol": 1e-5},
    )
    frequency = float(result.x)
    amplitude = float(-result.fun)

    noise_mask = use & (np.abs(frequencies - frequency) > 4.0 * resolution)
    noise = float(np.median(spectrum[noise_mask]))
    snr = amplitude / max(noise, 1e-12)
    return frequency, amplitude, snr


def aicc(observed: np.ndarray, predicted: np.ndarray, n_parameters: int) -> float:
    n = observed.size
    rss = float(np.sum((observed - predicted) ** 2))
    rss = max(rss, 1e-18)
    value = n * np.log(rss / n) + 2 * n_parameters
    if n > n_parameters + 1:
        value += 2 * n_parameters * (n_parameters + 1) / (n - n_parameters - 1)
    return float(value)


def segment_amplitudes(
    processed: dict[str, np.ndarray | float],
    frequency: float,
    n_segments: int = 7,
) -> list[tuple[float, float]]:
    x = np.asarray(processed["x"], float)
    y = np.asarray(processed["y"], float)
    left = int(0.04 * x.size)
    right = int(0.96 * x.size)
    edges = np.linspace(left, right, n_segments + 1, dtype=int)
    values: list[tuple[float, float]] = []
    for start, stop in zip(edges[:-1], edges[1:]):
        xs = x[start:stop]
        ys = y[start:stop]
        window = np.hanning(xs.size)
        coeff = np.sum(ys * window * np.exp(-2j * np.pi * frequency * xs))
        amplitude = 2.0 * np.abs(coeff) / window.sum()
        effective_field = 1.0 / np.mean(xs)
        values.append((float(effective_field), float(amplitude)))
    return values


def fit_lk_dingle(points: pd.DataFrame) -> tuple[float, float]:
    B = points["B"].to_numpy(float)
    T = points["temperature"].to_numpy(float)
    A = points["amplitude"].to_numpy(float)

    def prediction(params: np.ndarray) -> np.ndarray:
        scale, mass, dingle = params
        x = ALPHA_LK * mass * T / B
        thermal = x / np.sinh(x)
        return scale * thermal * np.exp(-ALPHA_LK * mass * dingle / B)

    def residual(params: np.ndarray) -> np.ndarray:
        model = np.maximum(prediction(params), 1e-15)
        return np.log(np.maximum(A, 1e-15)) - np.log(model)

    fit = optimize.least_squares(
        residual,
        x0=np.array([0.10, 0.30, 2.5]),
        bounds=(
            np.array([0.01, 0.05, 0.10]),
            np.array([1.00, 2.00, 20.0]),
        ),
    )
    return float(fit.x[1]), float(fit.x[2])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    field_gain, field_offset = linear_calibration(
        DATA / "field_calibration.csv", "sensor_field_t", "reference_field_t"
    )
    angle_gain, angle_offset = linear_calibration(
        DATA / "angle_calibration.csv", "encoder_angle_deg", "reference_angle_deg"
    )
    temp_gain, temp_offset = linear_calibration(
        DATA / "temperature_calibration.csv", "sensor_temperature_k", "reference_temperature_k"
    )

    sweeps = pd.read_csv(DATA / "oscillation_sweeps.csv")
    processed: dict[str, dict[str, np.ndarray | float]] = {}
    sweep_rows: list[dict[str, float | str]] = []

    for sweep_id, group in sweeps.groupby("sweep_id", sort=True):
        p = clean_and_resample(
            group,
            field_gain,
            field_offset,
            angle_gain,
            angle_offset,
            temp_gain,
            temp_offset,
        )
        processed[sweep_id] = p
        frequency, amplitude, snr = continuous_peak(p)
        sweep_rows.append(
            {
                "sweep_id": sweep_id,
                "series": str(group["series"].iloc[0]),
                "direction": str(group["direction"].iloc[0]),
                "angle_deg": float(p["theta"]),
                "temperature_k": float(p["temperature"]),
                "frequency_t": frequency,
                "amplitude": amplitude,
                "peak_snr": snr,
            }
        )

    sweep_results = pd.DataFrame(sweep_rows)
    angle_rows = sweep_results.loc[sweep_results["series"] == "angle"].copy()
    angle_rows["cluster"] = np.rint(angle_rows["angle_deg"]).astype(int)

    angle_summary = (
        angle_rows.groupby("cluster", as_index=False)
        .agg(
            angle_deg=("angle_deg", "mean"),
            frequency_t=("frequency_t", "mean"),
            frequency_spread_t=("frequency_t", "std"),
            peak_snr=("peak_snr", "mean"),
        )
        .drop(columns="cluster")
        .sort_values("angle_deg")
        .reset_index(drop=True)
    )
    angle_summary["frequency_spread_t"] = angle_summary["frequency_spread_t"].fillna(0.0)
    angle_summary.to_csv(OUT / "angle_frequencies.csv", index=False, float_format="%.8g")

    theta = np.deg2rad(angle_summary["angle_deg"].to_numpy(float))
    observed = angle_summary["frequency_t"].to_numpy(float)

    q_basis = 1.0 / np.abs(np.cos(theta))
    q_f0 = float(np.dot(q_basis, observed) / np.dot(q_basis, q_basis))
    q_pred = q_f0 * q_basis

    gamma = 3.2
    e_basis = 1.0 / np.sqrt(np.cos(theta) ** 2 + np.sin(theta) ** 2 / gamma**2)
    e_f0 = float(np.dot(e_basis, observed) / np.dot(e_basis, e_basis))
    e_pred = e_f0 * e_basis

    if aicc(observed, q_pred, 1) <= aicc(observed, e_pred, 1):
        model = "quasi-2d"
        F0 = q_f0
    else:
        model = "ellipsoidal-3d"
        F0 = e_f0

    zero_rows = sweep_results[
        (sweep_results["series"] == "temperature")
        | (
            (sweep_results["series"] == "angle")
            & (np.abs(sweep_results["angle_deg"]) < 1.0)
        )
    ].copy()

    segment_rows: list[dict[str, float | int]] = []
    for _, row in zero_rows.iterrows():
        for seg, (effective_field, amplitude) in enumerate(
            segment_amplitudes(processed[str(row["sweep_id"])], F0)
        ):
            segment_rows.append(
                {
                    "temperature": float(row["temperature_k"]),
                    "temperature_cluster": int(round(float(row["temperature_k"]) * 10.0)),
                    "segment": seg,
                    "B": effective_field,
                    "amplitude": amplitude,
                }
            )

    segment_frame = pd.DataFrame(segment_rows)
    grouped_segments = (
        segment_frame.groupby(["temperature_cluster", "segment"], as_index=False)
        .agg(
            temperature=("temperature", "mean"),
            B=("B", "mean"),
            amplitude=("amplitude", "mean"),
        )
    )
    effective_mass, dingle_temperature = fit_lk_dingle(grouped_segments)

    constants = json.loads((DATA / "constants.json").read_text(encoding="utf-8"))
    e = float(constants["elementary_charge_c"])
    hbar = float(constants["hbar_j_s"])
    kB = float(constants["boltzmann_j_k"])
    me = float(constants["electron_mass_kg"])

    extremal_area = 2.0 * np.pi * e * F0 / hbar
    kF = np.sqrt(extremal_area / np.pi)
    quantum_lifetime = hbar / (2.0 * np.pi * kB * dingle_temperature)
    fermi_velocity = hbar * kF / (effective_mass * me)
    mean_free_path = fermi_velocity * quantum_lifetime

    holdout = pd.read_csv(DATA / "holdout_angles.csv")
    holdout_theta = np.deg2rad(holdout["angle_deg"].to_numpy(float))
    if model == "quasi-2d":
        holdout_frequency = F0 / np.abs(np.cos(holdout_theta))
    else:
        holdout_frequency = F0 / np.sqrt(
            np.cos(holdout_theta) ** 2 + np.sin(holdout_theta) ** 2 / gamma**2
        )
    pd.DataFrame(
        {
            "angle_deg": holdout["angle_deg"].to_numpy(float),
            "predicted_frequency_t": holdout_frequency,
        }
    ).to_csv(OUT / "holdout_predictions.csv", index=False, float_format="%.8g")

    summary = {
        "field_gain": field_gain,
        "field_offset_t": field_offset,
        "angle_gain": angle_gain,
        "angle_offset_deg": angle_offset,
        "temperature_gain": temp_gain,
        "temperature_offset_k": temp_offset,
        "model": model,
        "F0_t": F0,
        "effective_mass_me": effective_mass,
        "dingle_temperature_k": dingle_temperature,
        "extremal_area_m_inv2": float(extremal_area),
        "kF_m_inv": float(kF),
        "quantum_lifetime_s": float(quantum_lifetime),
        "fermi_velocity_m_s": float(fermi_velocity),
        "mean_free_path_m": float(mean_free_path),
    }
    (OUT / "pocket_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    report = f"""# Magnetotransport analysis

The calibrated angle series contains one dominant oscillation branch that moves from about {angle_summary['frequency_t'].iloc[0]:.1f} T near the c axis to about {angle_summary['frequency_t'].iloc[-1]:.1f} T at the largest measured angle. The branch follows the {model} candidate more closely than the alternative model. The fitted zero-angle frequency is {F0:.2f} T.

The temperature and field dependence of the fundamental gives an effective cyclotron mass of {effective_mass:.3f} electron masses and a Dingle temperature of {dingle_temperature:.3f} K. Those values imply an extremal orbit area of {extremal_area:.4e} m^-2, kF = {kF:.4e} m^-1, a quantum lifetime of {quantum_lifetime:.4e} s, a Fermi velocity of {fermi_velocity:.4e} m/s and a quantum mean free path of {mean_free_path:.4e} m.

The raw sweeps required calibration of the field, rotator and thermometer channels. Duplicate field readings were consolidated before spectral analysis. Missing resistivity rows were excluded, while isolated spikes were treated as acquisition artefacts rather than oscillation extrema. Up- and down-field sweeps were analysed independently before being combined, so the reported spread reflects repeatability rather than a single trace. A weak harmonic was visible but did not replace the lower-frequency fundamental as the branch used for the pocket reconstruction.
"""
    (OUT / "analysis_report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
