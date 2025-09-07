
# Power Electronics – Laboratory Exercises

**Authors:**
* Andreas Karogiannis

## Overview

This repository contains the first four laboratory exercises for **Power Electronics**, focusing on signal analysis, rectifiers, DC-DC converters, and inverters.
The main goal is to understand the relationship between signal distortion, frequency content, duty cycle, firing angle, PI control, and power efficiency.

**Covered Topics:**

* Lab 1: Signal analysis using Fourier, THD, DPF, RMS
* Lab 2: Single-phase and three-phase rectifiers with RL loads and PI current control
* Lab 3: Buck converters and DC motor control via duty cycle and PI control
* Lab 4: Single-phase and three-phase inverters with PWM control

---

## Theoretical Background

* **Fourier & Signal Analysis:**

  * Discrete Fourier Transform (DFT) decomposes signals into frequency components.
  * THD, DPF, RMS, and power are calculated.
  * Higher distortion → higher THD, lower DPF & PF.

* **Rectifiers (Lab 2):**

  * Mono- and three-phase rectifiers feed RL loads.
  * Output current depends on load inductance L and firing angle α.
  * PI controller stabilizes current at a desired reference.

* **Buck Converters & DC Motors (Lab 3):**

  * Output voltage and current controlled via duty cycle.
  * PI controller ensures motor speed stability.
  * Higher L → smoother current and voltage.

* **Inverters (Lab 4):**

  * Single- and three-phase bridge inverters with square-wave outputs.
  * Firing angle α and PWM control harmonic content.
  * Higher switching frequency (mf) → denser waveforms, fewer harmonics.

---

## Simulations

* Implemented in **MATLAB/Python**.
* Parameters varied for each lab: L, R, α, duty cycle, TL (load torque), ma, mf.
* Graphs generated for voltage, current, harmonics, THD, DPF, RMS, and power.

---

## Key Results & Observations

1. **Signals & Fourier Analysis:**

   * More distorted signals → higher THD, lower DPF/PF.
   * RMS from spectrum ≈ RMS from definition.

2. **Rectifiers:**

   * Increasing L → output current smoother (more DC-like).
   * Increasing α → delayed conduction, lower current amplitude.

3. **Buck Converters & DC Motors:**

   * Higher duty cycle → higher voltage/current.
   * PI control stabilizes motor speed despite load changes.

4. **Inverters:**

   * Lower α → smoother waveforms, improved power factor.
   * Higher mf → more compact waveforms, fewer harmonics.

---

## Usage

1. Run the MATLAB/Python scripts for each lab.
2. Adjust parameters (L, R, α, duty cycle, TL, ma, mf) as desired.
3. Observe voltage, current, and harmonic waveforms.
4. Calculate THD, DPF, RMS, and power when required.

---


