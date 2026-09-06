Low-temperature magnetotransport measurements from a layered metal are in `/app/data/oscillation_sweeps.csv`. The field, rotator angle and thermometer channels each have separate calibration records in `/app/data/field_calibration.csv`, `/app/data/angle_calibration.csv` and `/app/data/temperature_calibration.csv`. The measurement manifest, candidate angular-frequency models, physical constants and requested holdout angles are also in `/app/data/`.

The oscillatory signal is much smaller than the smooth magnetoresistance background. The run contains repeated up- and down-field sweeps, isolated acquisition spikes, duplicated field readings, missing resistivity values, a weak harmonic, several field angles and a temperature series near the crystallographic c axis.

Determine the dominant quantum-oscillation branch, decide which candidate angular model is supported by the calibrated measurements, and extract the pocket parameters implied by the data.

Write `/app/output/angle_frequencies.csv` with columns `angle_deg`, `frequency_t`, `frequency_spread_t` and `peak_snr`. Write `/app/output/holdout_predictions.csv` with columns `angle_deg` and `predicted_frequency_t`.

Write `/app/output/pocket_summary.json` containing `field_gain`, `field_offset_t`, `angle_gain`, `angle_offset_deg`, `temperature_gain`, `temperature_offset_k`, `model`, `F0_t`, `effective_mass_me`, `dingle_temperature_k`, `extremal_area_m_inv2`, `kF_m_inv`, `quantum_lifetime_s`, `fermi_velocity_m_s` and `mean_free_path_m`.

Also write `/app/output/analysis_report.md`, briefly explaining the model decision, the parameter estimates and any data-quality choices that materially affected the result.
