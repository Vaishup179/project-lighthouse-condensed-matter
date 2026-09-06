# Solution explanation

## Research workflow represented by the task

This task mirrors a common condensed-matter transport workflow: taking raw magnetoresistance sweeps from a rotatable cryogenic measurement and turning the oscillatory component into a physically consistent description of a Fermi-surface pocket. The difficult part is not a single formula. The result depends on several decisions that interact with one another. The field, sample angle and temperature channels are all slightly miscalibrated; the magnetoresistance background is much larger than the quantum oscillation; the raw table contains duplicate rows, missing resistivity values and isolated spikes; a second harmonic is present; and the main physical parameters must agree across angle, temperature and field dependence.

The input data are synthetic but intentionally shaped like experimental data rather than a clean textbook signal. They are original to this task package, so no external measurement dataset or licence is required. The scientific relationships used to interpret the oscillations follow standard quantum-oscillation analysis.

## Interpreting and calibrating the measurements

An expert would begin by fitting the three calibration tables. Each table gives a measured instrument channel against a trusted reference value. A linear calibration is appropriate here because the residual errors are small and the reference points cover the operating range of the experiment. Those fits provide the gain and offset needed to convert the raw field sensor, rotator encoder and thermometer readings into calibrated physical values.

The raw sweep table should then be separated by sweep identifier. Missing resistivity values cannot contribute to the oscillatory analysis and should be removed. Exact duplicate field readings should be consolidated rather than treated as independent samples. Isolated spikes are acquisition artefacts; they should not be allowed to create large Fourier components. At the same time, aggressive smoothing would damage the oscillatory signal, so the cleaning has to be local and conservative.

The smooth magnetoresistance background is many orders larger than the oscillatory correction in relative terms. A low-order smooth background model is therefore removed independently from each sweep. The exact background estimator is not part of the grading contract; another scientifically reasonable smooth estimator can work as long as it preserves the oscillation frequency and damping. After background removal, the independent variable is changed from field to inverse field, because quantum oscillations are approximately periodic in \(1/B\). The data are then sampled on a uniform inverse-field grid before spectral estimation.

## Recovering the oscillation branch

Each calibrated angle has both up- and down-field traces. They are analysed separately first, which provides a useful repeatability check and prevents one bad trace from silently controlling the answer. A windowed spectrum in inverse field reveals a dominant fundamental peak and a weaker second harmonic. The lower-frequency branch is followed across angle. The harmonic should remain identifiable as a multiple of the same branch rather than being mistaken for a second Fermi-surface pocket.

The extracted fundamental frequency increases strongly with angle. Two candidate angular models are supplied with the task. For the quasi-two-dimensional model the frequency scales as \(F_0/|\cos	heta|\). The alternative is a three-dimensional ellipsoid with a fixed anisotropy supplied from an independent measurement. The calibrated frequencies follow the quasi-two-dimensional relation much more closely, so the fitted zero-angle frequency \(F_0\) is taken from that model. The reference solution obtains about 244 T.

This angular result is also used for the requested holdout predictions. Those predictions are not separate free parameters: they must follow from the selected model and the fitted \(F_0\).

## Effective mass and Dingle temperature

The temperature series near the crystallographic c axis contains information that is independent of the angular frequency fit. The fundamental amplitude decreases with temperature through the Lifshitz-Kosevich thermal reduction factor. It also grows toward higher magnetic field because of Dingle damping. Rather than using a single field value for the entire trace, the reference solution measures the fundamental amplitude in several inverse-field windows and fits the combined field-and-temperature dependence.

For a cyclotron mass \(m^\*\), temperature \(T\) and field \(B\), the thermal factor is

\[
R_T = X/\sinh X,
\]

with \(X\) proportional to \(m^\*T/B\). The field damping is represented by the standard Dingle factor. Fitting the amplitudes across several temperatures and field windows gives both the cyclotron effective mass and the Dingle temperature. The exact window layout is not graded; the verifier reconstructs the same physical quantities independently with a different segmentation. This is deliberate. A valid solution should survive a reasonable change in implementation details.

## Derived Fermi-surface quantities

Once \(F_0\), \(m^\*\) and the Dingle temperature are available, the remaining quantities are linked by standard physics and should be mutually consistent. Onsager's relation connects \(F_0\) to the extremal orbit area \(A_k\),

\[
A_k = 2\pi eF_0/\hbar.
\]

Assuming a circular in-plane orbit, \(k_F=\sqrt{A_k/\pi}\). The Dingle temperature gives the quantum lifetime through \(	au_q=\hbar/(2\pi k_B T_D)\). The Fermi velocity is \(v_F=\hbar k_F/(m^\*m_e)\), and the corresponding quantum mean free path is \(v_F	au_q\). These are downstream quantities, so an incorrect frequency, mass or Dingle temperature propagates into several later outputs. That cascading structure is part of the intended task difficulty.

## Verification strategy

The verifier does not compare the submission with a stored answer key. It loads the raw calibration and sweep files and reconstructs the scientific results independently. Calibration coefficients are refitted from the references. Quantum-oscillation frequencies are recovered from the raw sweeps using a separately written cleaning and inverse-field spectral pipeline. The angular model and \(F_0\) are determined from those independently recovered frequencies. Effective mass and Dingle temperature are estimated again from independently segmented amplitude data. The verifier then checks that the submitted derived quantities satisfy the physical relationships among \(F_0\), orbit area, \(k_F\), lifetime, Fermi velocity and mean free path.

The numerical tolerances are intentionally wider than the difference between reasonable signal-processing choices. This allows alternative valid background subtraction, interpolation and peak-refinement methods while still rejecting wrong harmonics, uncalibrated axes, incorrect model choices or hardcoded guesses. The holdout predictions are checked against the submitted model and against the raw-data reconstruction of \(F_0\).

The expected chain of work contains more than 25 meaningful decisions and computations: three calibrations, file and sweep inspection, missing-data handling, duplicate consolidation, spike treatment, background estimation, inverse-field conversion, resampling, spectral estimation, peak refinement, fundamental-versus-harmonic identification, direction-to-direction validation, angle aggregation, two model evaluations, model choice, zero-angle fit, temperature grouping, field-window amplitude extraction, coupled thermal/Dingle fitting, seven derived physical calculations and output/report generation. An early calibration or branch-selection error therefore damages several later quantities rather than remaining isolated.

## References

Onsager, L. (1952). Interpretation of the de Haas-van Alphen effect. *The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science, 43*(344), 1006-1008. https://doi.org/10.1080/14786440908521019  
[Provides the frequency-to-extremal-area relation used for the Fermi-surface reconstruction.]

Lifshitz, I. M., & Kosevich, A. M. (1956). Theory of magnetic susceptibility in metals at low temperature. *Soviet Physics JETP, 2*(4), 636-645.  
[Provides the thermal damping framework used for the cyclotron effective-mass inference.]

Shoenberg, D. (1984). *Magnetic Oscillations in Metals*. Cambridge University Press. https://doi.org/10.1017/CBO9780511897870  
[Provides the standard experimental and theoretical framework for quantum-oscillation spectra, Dingle damping and Fermi-surface interpretation.]
