# Scenario: Earth at Venus's Orbit with a 90-Bar Dry-Air Atmosphere (HITRAN-Supported)

This directory contains the physical model, calculations, and simulation code for a hypothetical scenario where Earth is moved to Venus's orbit (0.723 AU) and its atmospheric mass is increased to reach a surface pressure ($P_s$) of 90 bar. The atmospheric composition is kept strictly identical to modern dry air ($\sim 78\% \text{ N}_2$, $\sim 21\% \text{ O}_2$, $\sim 1\% \text{ Ar}$), and the oceans are absent.

The calculations are supported directly by the **HITRAN database** for collision-induced absorption (CIA).

## The Original Prompt / Scenario

> **Planetary Scenario:** Earth is moved to Venus's orbit.
> 
> ### Core Astronomical & Orbital Parameters
> * **Planetary Distance:** $d = 0.723 \text{ AU}$
> * **Base Solar Constant at 1 AU:** $S_0 = 1361 \text{ W/m}^2$
> * **Gravitational Acceleration:** $g = 9.81 \text{ m/s}^2$
> 
> ### The Scenario: The 90-Bar Dry-Air Atmosphere
> Earth's atmospheric mass is increased to reach a surface pressure ($P_s$) of 90 bar at sea level, but its composition is kept strictly identical to modern dry air ($\sim 78\% \text{ N}_2$, $\sim 21\% \text{ O}_2$, $\sim 1\% \text{ Ar}$). Earth's oceans are completely suppressed or absent.
> 
> ### Physical Mechanisms:
> 1. **Solar Forcing & Albedo:** Calculate the incident solar radiation at $0.723 \text{ AU}$. Account for the extreme Rayleigh scattering caused by a 90-bar diatomic column, which significantly raises the planetary spherical albedo ($\alpha$) to a range of $0.6 \text{ to } 0.7$. Use $\alpha = 0.7$ for consistency. Calculate the resulting planetary equilibrium temperature ($T_{eq}$).
> 2. **Infrared Opacity via Collision-Induced Absorption (CIA):** While $\text{N}_2$ and $\text{O}_2$ are transparent to thermal IR at 1 bar, explicitly account for the fact that at 90 bar, the probability of collisions induces temporary dipole moments. Because CIA opacity scales non-linearly with the square of the density ($\rho^2$), the atmosphere becomes opaque to thermal IR. Identify the effective emission level and surface temperature.

---

## Physical Derivation & Calculations

### 1. Solar Forcing & Equilibrium Temperature ($T_{eq}$)

The solar flux $S$ at a distance $d = 0.723 \text{ AU}$ is determined by the inverse-square law:
$$S = S_0 \left(\frac{1\text{ AU}}{d}\right)^2 = 1361 \text{ W/m}^2 \times \left(\frac{1}{0.723}\right)^2 \approx 2603.64 \text{ W/m}^2$$

With a planetary spherical albedo $\alpha = 0.7$, the planet reflects $70\%$ of the incoming radiation back to space, primarily due to intense Rayleigh scattering by the dense 90-bar diatomic atmosphere. The net absorbed solar radiation per unit surface area (averaged over the sphere) is:
$$F_{\text{abs}} = \frac{S (1 - \alpha)}{4} = \frac{2603.64 \times (1 - 0.7)}{4} \approx 195.27 \text{ W/m}^2$$

At thermal equilibrium, the outgoing longwave radiation (OLR) must balance the net absorbed solar radiation:
$$\text{OLR} = F_{\text{abs}} = 195.27 \text{ W/m}^2$$

The planetary equilibrium temperature ($T_{eq}$), representing the effective temperature at which the planet radiates as a blackbody to space, is:
$$T_{eq} = \left(\frac{\text{OLR}}{\sigma}\right)^{1/4} = \left(\frac{195.27 \text{ W/m}^2}{5.6704 \times 10^{-8} \text{ W/m}^2/\text{K}^4}\right)^{1/4} \approx 242.25 \text{ K} \quad (-30.90\text{ }^\circ\text{C})$$

---

### 2. Vertical Atmospheric Structure

To resolve the temperature profile, we construct a 1D radiative-convective equilibrium model. The atmosphere is divided into two primary zones:
1. **The Stratosphere (Radiative Equilibrium):** Above the tropopause ($P \le P_t$), heat transport is purely radiative. Under the grey-atmosphere Eddington approximation, the temperature profile is:
   $$T^4(\tau) = T_{eq}^4 \left(\frac{3}{4}\tau + \frac{1}{2}\right)$$
   where $\tau$ is the infrared optical depth measured downward from the top of the atmosphere.
2. **The Troposphere (Convective Equilibrium):** Below the tropopause ($P > P_t$), the steep radiative temperature gradient is unstable to convection. Because the atmosphere is dry (no water vapor or condensation), the temperature profile follows the dry adiabat:
   $$T(P) = T_t \left(\frac{P}{P_t}\right)^\kappa$$
   where $\kappa = \frac{R_d}{C_p}$ is the adiabatic index, $T_t$ is the tropopause temperature, and $P_t$ is the tropopause pressure.

#### Atmospheric Properties of Dry Air:
* **Mean Molecular Weight ($M$):** For $\sim 78\%\text{ N}_2$, $\sim 21\%\text{ O}_2$, $\sim 1\%\text{ Ar}$, and trace $\text{CO}_2$:
  $$M \approx 28.965 \text{ g/mol} = 0.028965 \text{ kg/mol}$$
* **Specific Gas Constant ($R_d$):**
  $$R_d = \frac{R}{M} \approx 287.05 \text{ J/(kg K)}$$
* **Specific Heat Capacity ($C_p$):** For diatomic molecules $\text{N}_2$ and $\text{O}_2$, the molar heat capacity is $C_{p,m} = \frac{7}{2}R$; for monatomic $\text{Ar}$, it is $C_{p,m} = \frac{5}{2}R$.
  $$C_p \approx 1002.05 \text{ J/(kg K)}$$
* **Adiabatic Index ($\kappa$):**
  $$\kappa = \frac{R_d}{C_p} \approx 0.2865$$

---

### 3. Spectroscopic Opacity from the HITRAN Database

Because $\text{N}_2$ and $\text{O}_2$ are homonuclear diatomic molecules, they lack permanent dipole moments and are transparent to infrared radiation at low pressures. However, in a dense 90-bar atmosphere, molecular collisions are frequent enough that the transient distortion of their electron clouds induces temporary dipole moments. This gives rise to **Collision-Induced Absorption (CIA)**, which is proportional to the square of the density ($\rho^2$).

The mass absorption coefficient $\kappa_{\text{IR}}$ is:
$$\kappa_{\text{IR}}(T, \rho) = \kappa_{R,0}(T) \cdot \rho$$
where $\kappa_{R,0}(T)$ is the temperature-dependent density-independent mean absorption coefficient (in $\text{m}^5 \text{ kg}^{-2}$).
Using the ideal gas law ($\rho = \frac{P}{R_d T}$), the differential optical depth $d\tau$ is:
$$d\tau = \frac{\kappa_{\text{IR}}}{g} dP = \frac{\kappa_{R,0}(T) P}{g R_d T} dP$$

Integrating from the top of the atmosphere ($P=0$) down to a pressure level $P$:
$$\tau(P) = \int_0^P \frac{\kappa_{R,0}(T(P')) P'}{g R_d T(P')} dP'$$

We retrieve the raw collision-induced absorption (CIA) data files from the **HITRAN database** for $\text{N}_2$-$\text{N}_2$ (`N2-N2_2021.cia`), $\text{O}_2$-$\text{N}_2$ (`O2-N2_2024.cia`), and $\text{O}_2$-$\text{O}_2$ (`O2-O2_2024.cia`). We group the data blocks into separate bands and perform temperature-dependent linear interpolation for each band. The combined dry-air CIA spectrum is calculated as:
$$\text{CIA}_{\text{air}}(\nu, T) = f_{N_2}^2 \text{CIA}_{N_2-N_2}(\nu, T) + f_{O_2} f_{N_2} \text{CIA}_{O_2-N_2}(\nu, T) + f_{O_2}^2 \text{CIA}_{O_2-O_2}(\nu, T)$$

We calculate two types of mean opacities:
1. **Planck Mean (Case 1):** The arithmetic mean weighted by the Planck function, representing the total integrated energy trapping capability:
   $$\kappa_{\text{Planck}}(T) = \left(\frac{N_A}{M}\right)^2 \int_0^\infty \text{CIA}_{\text{air}}(\nu, T) \frac{15}{\pi^4} \frac{x^3}{e^x - 1} \frac{hc}{k_B T} d\nu$$
2. **Rosseland Mean (Case 2):** The harmonic mean weighted by the temperature derivative of the Planck function, which represents the radiative diffusion through the optically thick core. To avoid numerical singularities at the band edges, we integrate over the active absorbing far-infrared band ($0 - 650\text{ cm}^{-1}$):
   $$\frac{1}{\kappa_{\text{Rosseland}}(T)} = \left(\frac{M}{N_A}\right)^2 \int_{0}^{650} \frac{w_x(x)}{\text{CIA}_{\text{air}}(\nu, T)} d\nu$$
   where $w_x(x)$ is the normalized Rosseland weight.

---

## 4. Quantitative Results from Numerical Model

Solving this radiative-convective system iteratively yields the following parameters:

| Parameter | Planck Mean (HITRAN) | Rosseland Mean (HITRAN) | Baseline Constant ($10^{-4}$) |
| :--- | :--- | :--- | :--- |
| **Surface Temperature ($T_s$)** | **$720.1\text{ K}$ ($447.0\text{ }^\circ\text{C}$)** | **$511.7\text{ K}$ ($238.5\text{ }^\circ\text{C}$)** | **$830.4\text{ K}$ ($557.3\text{ }^\circ\text{C}$)** |
| **Total Optical Depth ($\tau_s$)** | $402.2$ | $62.4$ | $1210.4$ |
| **Effective Emission Pressure ($P_e$)** | $1.802\text{ bar}$ | $6.439\text{ bar}$ | $1.172\text{ bar}$ |
| **Effective Emission Height ($z_e$)** | $49.63\text{ km}$ | $27.72\text{ km}$ | $60.37\text{ km}$ |
| **Effective Emission Temp ($T_e$)** | $242.8\text{ K}$ ($-30.4\text{ }^\circ\text{C}$) | $242.4\text{ K}$ ($-30.8\text{ }^\circ\text{C}$) | $241.5\text{ K}$ ($-31.7\text{ }^\circ\text{C}$) |
| **Tropopause Pressure ($P_t$)** | $4.323\text{ bar}$ | $8.735\text{ bar}$ | $1.622\text{ bar}$ |
| **Tropopause Height ($z_t$)** | $42.73\text{ km}$ | $25.47\text{ km}$ | $57.98\text{ km}$ |
| **Tropopause Temperature ($T_t$)** | $301.8\text{ K}$ ($28.7\text{ }^\circ\text{C}$) | $262.3\text{ K}$ ($-10.9\text{ }^\circ\text{C}$) | $262.8\text{ K}$ ($-10.4\text{ }^\circ\text{C}$) |

> The effective emission temperature $T_e$ recovers $T_{eq} \approx 242.25\text{ K}$ in every case. This is a built-in consistency check: in the grey Eddington model $T^4(\tau{=}2/3) = T_{eq}^4(\tfrac34\cdot\tfrac23 + \tfrac12) = T_{eq}^4$, so the emission level *must* sit at $T_{eq}$ regardless of the opacity used. It validates the optical-depth integration, not the surface temperature.

### Interpretation

The choice of mean opacity is not a free parameter — it determines which physical regime the model represents, and the two means bracket the answer for very different reasons.

* **Rosseland Mean — the physically appropriate estimate for the deep atmosphere.** The Rosseland mean is the harmonic average that governs radiative *diffusion* through an optically thick medium, which is exactly the regime of the deep atmosphere here ($\tau_s \gg 1$). Because it is a harmonic average, it is dominated by the most transparent parts of the spectrum: thermal radiation preferentially leaks out through the spectral *windows*, so those windows set the true energy loss rate. This yields a smaller effective optical depth ($\tau_s \approx 62$), an emission level deep in the atmosphere ($P_e \approx 6.4\text{ bar}$, $\approx 28\text{ km}$), and a surface temperature of **$\approx 512\text{ K}$**. For a grey-diffusion idealization of this scenario, this is the estimate to trust.

* **Planck Mean — an upper bound that over-states trapping.** The Planck mean is an arithmetic average weighted by the emitted spectrum; it is the correct mean only in the optically *thin* limit. Applied to a strongly non-grey, optically thick atmosphere it over-weights the strong absorption bands and ignores the fact that radiation escapes through the windows, inflating the effective optical depth to $\tau_s \approx 402$ and the surface temperature to **$\approx 720\text{ K}$**. The numerical closeness of this figure to Venus's observed $\approx 740\text{ K}$ is **coincidental** and should not be read as validation: real Venus is heated by a massive $\text{CO}_2$ atmosphere with $\text{SO}_2$/cloud opacity closing the windows, a mechanism absent from this dry $\text{N}_2$/$\text{O}_2$ model. Treat $720\text{ K}$ as a loose upper bound, not a prediction.

* **Bottom line.** A 90-bar dry-air atmosphere at Venus's orbit produces a strong CIA greenhouse — surface temperatures of several hundred Kelvin above $T_{eq}$ — but the defensible range is roughly **$510\text{–}600\text{ K}$** (Rosseland to intermediate), not a clean Venus match. The headline result is that N₂/O₂ collision-induced absorption *alone* can lift the surface well above the boiling point of water, not that it reproduces Venus.

> **Modeling caveats.** This is a grey, single-mean-opacity radiative-convective model: neither the Planck nor the Rosseland mean is exactly correct for a non-grey spectrum, and a true band-by-band (or correlated-$k$) calculation would fall between them. The infrared optical depth is integrated with a cumulative trapezoid rule from $P=0$; solar absorption within the atmosphere, scattering-greenhouse coupling, and lapse-rate feedbacks are not modeled.

---

## 5. Spectral Coverage & Limitations

The radiative transfer here is **not** frequency-resolved. The full infrared spectrum is collapsed into a *single* mean opacity per temperature (Planck or Rosseland) before the grey two-stream relation $T^4(\tau) = T_{eq}^4(\tfrac34\tau + \tfrac12)$ is applied. Frequencies enter only when that mean is computed over the grid `nu_grid = 1–3000 cm⁻¹`. Readers should be aware of exactly which absorption this captures — and what it leaves out.

### What opacity is included
**Only $\text{N}_2$/$\text{O}_2$ collision-induced absorption (CIA)** from the three HITRAN files. The bands that fall in the thermal range actually used by the model are:

| Band | Range (cm⁻¹) | In model? |
| :--- | :--- | :--- |
| $\text{N}_2$–$\text{N}_2$ rototranslational | $0$–$450$ | ✓ |
| $\text{O}_2$–$\text{O}_2$ / $\text{O}_2$–$\text{N}_2$ collision | $1150$–$1950$ | ✓ |
| collision fundamental region | $1850$–$3000$ | ✓ |
| **(no absorber)** | **$450$–$1150$** | **✗ — zero opacity** |
| $\text{N}_2$–$\text{N}_2$ weak band | $4300$–$5000$ | ✗ (above the $3000\text{ cm}^{-1}$ grid cutoff) |
| $\text{O}_2$ near-IR / visible (Chappuis, Wulf, …) | $7500$–$33000$ | ✗ (not thermal; solar absorption is represented only through the albedo $\alpha$) |

### What is left out, and why it matters
* **A fully transparent window at $450$–$1150\text{ cm}^{-1}$ sits directly on the thermal peak.** This region has *no* absorber in the model, yet it carries roughly **$60\%$ of the blackbody emission at $240$–$300\text{ K}$** and still $\sim 20\%$ at $720\text{ K}$. It is the dominant escape route for outgoing longwave radiation, which is precisely why the harmonic (Rosseland) mean — physically appropriate for the optically thick interior — yields a much lower optical depth and surface temperature than the Planck mean.
* **$\text{CO}_2$-specific bands are entirely absent.** $\text{CO}_2$ enters the code only through the mean molecular weight $M$ and heat capacity $C_p$ (its $0.00036$ mole fraction); it contributes **zero opacity**. Its strong $15\,\mu\text{m}$ ($667\text{ cm}^{-1}$) bending band — which on real Earth and Venus partly closes exactly the window above — is not modeled. The same applies to $\text{H}_2\text{O}$ and $\text{O}_3$ bands. For this strictly dry $\text{N}_2$/$\text{O}_2$ composition that omission is *consistent* (there is no $\text{CO}_2$ band opacity to add at $360\text{ ppm}$ without line data), but it means the model cannot close the thermal window the way a real $\text{CO}_2$-rich atmosphere does.
* **The $3000\text{ cm}^{-1}$ grid cutoff** drops the weak $\text{N}_2$–$\text{N}_2$ $4300$–$5000\text{ cm}^{-1}$ band and up to $\sim 14\%$ of the Planck weight at $720\text{ K}$. Because that excluded region is essentially window (negligible CIA), the effect on the opacity is small, though it slightly biases the Planck-mean normalization at the highest temperatures.
* **The Rosseland mean is restricted to $\le 650\text{ cm}^{-1}$** (to avoid singularities in transparent windows), so it effectively samples *only* the far-IR $\text{N}_2$–$\text{N}_2$ rototranslational band — which holds $\sim 57\%$ of the emission at $242\text{ K}$ but only $\sim 7\%$ at $720\text{ K}$ — and ignores the $1150$–$3000\text{ cm}^{-1}$ collision bands.

**In short:** the model faithfully includes the non-$\text{CO}_2$ (N₂/O₂ collision) infrared bands available in HITRAN below $3000\text{ cm}^{-1}$, but it is a grey-mean approximation with a wide-open transparent window at the thermal peak and no $\text{CO}_2$/$\text{H}_2\text{O}$/$\text{O}_3$ line opacity. A frequency-resolved (band-by-band or correlated-$k$) treatment that included those line bands would be required to close the window and is the natural next step.

---

## Simulation Code

The code used to run these calculations and generate the temperature profile plot is in [climate_model.py](file:///home/cmerk/repos/maps/earth_on_venus/climate_model.py).

To execute the code and generate the plot:
```bash
uv run python earth_on_venus/climate_model.py
```

### Temperature Profile Plot
The model produces the following temperature profile as a function of pressure and altitude:

![Temperature Profile](./temperature_profile.png)
