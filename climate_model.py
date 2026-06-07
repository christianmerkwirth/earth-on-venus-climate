import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------------------------------------------------
# 1. HITRAN CIA PARSER AND INTERPOLATOR
# -------------------------------------------------------------
def parse_cia_file_by_bands(filepath):
    """
    Parses a HITRAN .cia file and groups the blocks into separate bands
    based on their wavenumber range (WMIN, WMAX).
    Returns a list of dicts: [band1_dict, band2_dict, ...]
    where each band_dict maps temperature -> (wavenumbers, cia_values)
    """
    bands = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    i = 0
    while i < len(lines):
        line = lines[i]
        if len(line.strip()) == 0 or line.startswith('!'):
            i += 1
            continue
        try:
            molec = line[0:20].strip()
            wmin = float(line[20:30])
            wmax = float(line[30:40])
            npt = int(line[40:47])
            tem = float(line[47:54])
            
            wavenumbers = []
            cia_values = []
            for _ in range(npt):
                i += 1
                parts = lines[i].split()
                wavenumbers.append(float(parts[0]))
                cia_values.append(float(parts[1]))
            
            wavenumbers = np.array(wavenumbers)
            cia_values = np.array(cia_values)
            
            # Group into bands
            found_band = False
            for band_dict in bands:
                ref_temp = list(band_dict.keys())[0]
                ref_wmin, ref_wmax = band_dict[ref_temp][2], band_dict[ref_temp][3]
                if abs(wmin - ref_wmin) < 50.0 and abs(wmax - ref_wmax) < 50.0:
                    band_dict[tem] = (wavenumbers, cia_values, wmin, wmax)
                    found_band = True
                    break
            
            if not found_band:
                new_band_dict = {tem: (wavenumbers, cia_values, wmin, wmax)}
                bands.append(new_band_dict)
                
            i += 1
        except Exception as e:
            i += 1
            continue
            
    return bands

def interpolate_band_to_temp(band_dict, T, nu_grid):
    """
    Interpolates a band's spectra to a target temperature T, and places it on nu_grid.
    """
    avail_temps = sorted(band_dict.keys())
    T_clipped = np.clip(T, avail_temps[0], avail_temps[-1])
    
    if len(avail_temps) == 1:
        nus, cia, _, _ = band_dict[avail_temps[0]]
        return np.interp(nu_grid, nus, cia, left=0.0, right=0.0)
        
    idx = np.searchsorted(avail_temps, T_clipped)
    if idx == 0:
        T1, T2 = avail_temps[0], avail_temps[0]
        w1, w2 = 0.5, 0.5
    elif idx == len(avail_temps):
        T1, T2 = avail_temps[-1], avail_temps[-1]
        w1, w2 = 0.5, 0.5
    else:
        T1, T2 = avail_temps[idx-1], avail_temps[idx]
        w2 = (T_clipped - T1) / (T2 - T1)
        w1 = 1.0 - w2
        
    nus1, cia1, _, _ = band_dict[T1]
    nus2, cia2, _, _ = band_dict[T2]
    
    cia1_interp = np.interp(nu_grid, nus1, cia1, left=0.0, right=0.0)
    cia2_interp = np.interp(nu_grid, nus2, cia2, left=0.0, right=0.0)
    
    return w1 * cia1_interp + w2 * cia2_interp

# -------------------------------------------------------------
# 2. RADIATIVE-CONVECTIVE EQUILIBRIUM SOLVER
# -------------------------------------------------------------
def solve_atmosphere(n2_n2_bands, o2_n2_bands, o2_o2_bands, mode="planck"):
    """
    Solves for the radiative-convective equilibrium temperature profile.
    mode can be: "planck" (Planck mean), "rosseland" (Rosseland mean), or "constant" (1e-4)
    """
    # Physical Constants
    S0 = 1361.0        # Solar constant at 1 AU (W/m^2)
    d = 0.723          # Planetary distance in AU (Venusian orbit)
    alpha = 0.7        # Spherical albedo
    sigma = 5.670374419e-8 # Stefan-Boltzmann constant (W/m^2/K^4)
    g = 9.81           # Gravitational acceleration (m/s^2)
    P_s = 90.0e5       # Surface pressure in Pa (90 bar)

    # Dry Air properties
    f_N2 = 0.78084
    f_O2 = 0.20946
    f_Ar = 0.00934
    f_CO2 = 0.00036

    M_N2, M_O2, M_Ar, M_CO2 = 0.028013, 0.031999, 0.039948, 0.04401
    M = f_N2*M_N2 + f_O2*M_O2 + f_Ar*M_Ar + f_CO2*M_CO2
    R_d = 8.314462618 / M
    
    Cp_m = (f_N2 + f_O2)*3.5*8.314462618 + f_Ar*2.5*8.314462618 + f_CO2*4.0*8.314462618
    Cp = Cp_m / M
    kappa = R_d / Cp

    # Solar forcing
    S = S0 / (d**2)
    F_abs = S * (1.0 - alpha) / 4.0
    T_eq = (F_abs / sigma)**0.25

    # Grid setup
    N = 2000
    P_grid = np.linspace(0.001e5, P_s, N)
    dP = P_grid[1] - P_grid[0]
    
    # Unified wavenumber grid for computing spectroscopic means
    nu_grid = np.linspace(1.0, 3000.0, 1500) # cm^-1
    
    # Precompute or define function to get k_mean (m^-1 amagat^-2) at a given T
    memoized_k = {}
    def get_k_mean(T):
        T_key = round(T, 1)
        if T_key in memoized_k:
            return memoized_k[T_key]
            
        if mode == "constant":
            return 1.0e-4
            
        # Sum up bands for all three pairs at temperature T
        cia_n2_n2 = np.zeros_like(nu_grid)
        for band in n2_n2_bands:
            cia_n2_n2 += interpolate_band_to_temp(band, T, nu_grid)
            
        cia_o2_n2 = np.zeros_like(nu_grid)
        for band in o2_n2_bands:
            cia_o2_n2 += interpolate_band_to_temp(band, T, nu_grid)
            
        cia_o2_o2 = np.zeros_like(nu_grid)
        for band in o2_o2_bands:
            cia_o2_o2 += interpolate_band_to_temp(band, T, nu_grid)
            
        cia_air = (f_N2**2)*cia_n2_n2 + f_O2*f_N2*cia_o2_n2 + (f_O2**2)*cia_o2_o2
        
        # Planck weight and Rosseland weight
        c2 = 1.438776877
        x = c2 * nu_grid / T
        
        if mode == "planck":
            B_w = (15.0 / np.pi**4) * (x**3) / (np.exp(x) - 1.0)
            planck_mean_cia = np.sum(0.5 * (cia_air[:-1]*B_w[:-1] + cia_air[1:]*B_w[1:]) * np.diff(x))
            k_val = planck_mean_cia * 7.218786e40
        else: # rosseland
            w_x = (15.0 / (4.0 * np.pi**4)) * (x**4) * np.exp(x) / (np.exp(x) - 1.0)**2
            
            # We integrate over the main absorbing region to avoid edge divisions by 0
            mask_band = (nu_grid <= 650.0)
            x_band = x[mask_band]
            cia_band = cia_air[mask_band]
            w_band = w_x[mask_band]
            
            mask_active = (cia_band > 1e-47)
            if np.any(mask_active):
                x_act = x_band[mask_active]
                cia_act = cia_band[mask_active]
                w_act = w_band[mask_active]
                w_act_norm = w_act / np.sum(0.5 * (w_act[:-1] + w_act[1:]) * np.diff(x_act))
                inv_rosseland = np.sum(0.5 * (w_act_norm[:-1]/cia_act[:-1] + w_act_norm[1:]/cia_act[1:]) * np.diff(x_act))
                rosseland_mean_cia = 1.0 / inv_rosseland
            else:
                rosseland_mean_cia = 1e-45
            k_val = rosseland_mean_cia * 7.218786e40
            
        memoized_k[T_key] = k_val
        return k_val

    # Initialize T(P) profile
    T_profile = np.ones(N) * T_eq
    
    # STP parameters for density-to-amagat conversion
    T_STP = 273.15
    P_STP = 101325.0
    rho_STP = P_STP / (R_d * T_STP)
    
    max_iter = 100
    tol = 1e-4
    for iteration in range(max_iter):
        T_old = T_profile.copy()
        
        # Calculate local k_mean for each grid cell and compute optical depth:
        # dtau = k_mean(T) * (P / (g * rho_STP^2 * R_d * T)) dP
        k_grid = np.array([get_k_mean(t) for t in T_profile])
        
        integrand = (k_grid * P_grid) / (g * (rho_STP**2) * R_d * T_profile)
        # Cumulative trapezoidal integration of optical depth. The integrand
        # vanishes at P=0, so the contribution from 0 -> P_grid[0] is the
        # triangle 0.5 * integrand[0] * P_grid[0]; subsequent intervals use the
        # trapezoid rule. This is more accurate than a forward-rectangle sum and
        # avoids ignoring the opacity above the first grid level.
        seg = 0.5 * (integrand[1:] + integrand[:-1]) * dP
        tau_grid = np.empty(N)
        tau_grid[0] = 0.5 * integrand[0] * P_grid[0]
        tau_grid[1:] = tau_grid[0] + np.cumsum(seg)
        
        # Radiative equilibrium profile
        T_rad = T_eq * (0.75 * tau_grid + 0.5)**0.25
        
        # Convective check: dT/dP = kappa * T / P
        dT_rad_dP = np.gradient(T_rad, P_grid)
        dT_ad_dP = kappa * T_rad / P_grid
        
        convective_mask = dT_rad_dP > dT_ad_dP
        if np.any(convective_mask):
            tropopause_idx = np.where(convective_mask)[0][0]
        else:
            tropopause_idx = N - 1
            
        P_t = P_grid[tropopause_idx]
        T_t = T_rad[tropopause_idx]
        
        T_profile[:tropopause_idx] = T_rad[:tropopause_idx]
        T_profile[tropopause_idx:] = T_t * (P_grid[tropopause_idx:] / P_t)**kappa
        
        diff = np.max(np.abs(T_profile - T_old))
        if diff < tol:
            break
            
    # Calculate geopotential heights
    z_grid = np.zeros(N)
    for i in range(N-2, -1, -1):
        P_avg = 0.5 * (P_grid[i] + P_grid[i+1])
        T_avg = 0.5 * (T_profile[i] + T_profile[i+1])
        dz = (R_d * T_avg / (g * P_avg)) * (P_grid[i+1] - P_grid[i])
        z_grid[i] = z_grid[i+1] + dz
        
    idx_emission = np.argmin(np.abs(tau_grid - 2.0/3.0))
    
    return {
        "P": P_grid,
        "T": T_profile,
        "z": z_grid,
        "tau": tau_grid,
        "P_t": P_t,
        "T_t": T_t,
        "z_t": z_grid[tropopause_idx],
        "P_e": P_grid[idx_emission],
        "T_e": T_profile[idx_emission],
        "z_e": z_grid[idx_emission],
        "T_s": T_profile[-1],
        "tau_s": tau_grid[-1]
    }

# -------------------------------------------------------------
# 3. MAIN RUN AND PLOTTING
# -------------------------------------------------------------
def run_model():
    # Make sure we run relative to the script's directory so paths are consistent
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dir_path = os.path.join(script_dir, "cia_data")
    os.makedirs(dir_path, exist_ok=True)
    
    # Download files if missing
    urls = {
        "N2-N2_2021.cia": "https://hitran.org/data/CIA/main/N2-N2_2021.cia",
        "O2-N2_2024.cia": "https://hitran.org/data/CIA/main/O2-N2_2024.cia",
        "O2-O2_2024.cia": "https://hitran.org/data/CIA/main/O2-O2_2024.cia"
    }
    import urllib.request
    for fname, url in urls.items():
        fpath = os.path.join(dir_path, fname)
        if not os.path.exists(fpath):
            print(f"Downloading {fname} from HITRAN...")
            urllib.request.urlretrieve(url, fpath)
            print(f"Downloaded {fname}")
            
    n2_n2_bands = parse_cia_file_by_bands(os.path.join(dir_path, "N2-N2_2021.cia"))
    o2_n2_bands = parse_cia_file_by_bands(os.path.join(dir_path, "O2-N2_2024.cia"))
    o2_o2_bands = parse_cia_file_by_bands(os.path.join(dir_path, "O2-O2_2024.cia"))
    
    print("Solving Case 1: Planck Mean (HITRAN-derived)...")
    res_planck = solve_atmosphere(n2_n2_bands, o2_n2_bands, o2_o2_bands, mode="planck")
    
    print("Solving Case 2: Rosseland Mean (HITRAN-derived)...")
    res_rosseland = solve_atmosphere(n2_n2_bands, o2_n2_bands, o2_o2_bands, mode="rosseland")
    
    print("Solving Case 3: Constant k_mean (1e-4 baseline)...")
    res_const = solve_atmosphere(n2_n2_bands, o2_n2_bands, o2_o2_bands, mode="constant")
    
    # Output comparison table
    print("\n" + "="*80)
    print("PLANETARY STRUCTURE COMPARISON (HITRAN SUPPORTED)")
    print("="*80)
    print(f"{'Parameter':<35} | {'Planck Mean':<12} | {'Rosseland':<12} | {'Baseline (1e-4)':<15}")
    print("-" * 83)
    print(f"{'Surface Temperature (T_s)':<35} | {res_planck['T_s']:10.1f} K | {res_rosseland['T_s']:10.1f} K | {res_const['T_s']:13.1f} K")
    print(f"{'':<35} | ({res_planck['T_s']-273.15:6.1f} °C) | ({res_rosseland['T_s']-273.15:6.1f} °C) | ({res_const['T_s']-273.15:9.1f} °C)")
    print(f"{'Total Optical Depth (tau_s)':<35} | {res_planck['tau_s']:10.1f}   | {res_rosseland['tau_s']:10.1f}   | {res_const['tau_s']:13.1f}  ")
    print(f"{'Effective Emission Pressure (P_e)':<35} | {res_planck['P_e']/1e5:10.3f} bar | {res_rosseland['P_e']/1e5:10.3f} bar | {res_const['P_e']/1e5:13.3f} bar")
    print(f"{'Effective Emission Height (z_e)':<35} | {res_planck['z_e']/1000:10.2f} km  | {res_rosseland['z_e']/1000:10.2f} km  | {res_const['z_e']/1000:13.2f} km ")
    print(f"{'Effective Emission Temp (T_e)':<35} | {res_planck['T_e']:10.1f} K  | {res_rosseland['T_e']:10.1f} K  | {res_const['T_e']:13.1f} K ")
    print(f"{'Tropopause Pressure (P_t)':<35} | {res_planck['P_t']/1e5:10.3f} bar | {res_rosseland['P_t']/1e5:10.3f} bar | {res_const['P_t']/1e5:13.3f} bar")
    print(f"{'Tropopause Height (z_t)':<35} | {res_planck['z_t']/1000:10.2f} km  | {res_rosseland['z_t']/1000:10.2f} km  | {res_const['z_t']/1000:13.2f} km ")
    print(f"{'Tropopause Temperature (T_t)':<35} | {res_planck['T_t']:10.1f} K  | {res_rosseland['T_t']:10.1f} K  | {res_const['T_t']:13.1f} K ")
    
    # -------------------------------------------------------------
    # PLOTTING COMPARISON
    # -------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Log Pressure vs Temperature
    ax1.plot(res_planck['T'], res_planck['P'] / 1e5, 'r-', lw=2, label="Planck Mean (HITRAN)")
    ax1.plot(res_rosseland['T'], res_rosseland['P'] / 1e5, 'b-', lw=2, label="Rosseland Mean (HITRAN)")
    ax1.plot(res_const['T'], res_const['P'] / 1e5, 'k--', lw=1.5, label="Constant 1e-4 Baseline")
    
    ax1.axhline(res_planck['P_e']/1e5, color='r', linestyle=':', alpha=0.7)
    ax1.axhline(res_rosseland['P_e']/1e5, color='b', linestyle=':', alpha=0.7)
    
    ax1.set_yscale('log')
    ax1.invert_yaxis()
    ax1.set_xlabel("Temperature (K)", fontsize=12)
    ax1.set_ylabel("Pressure (bar)", fontsize=12)
    ax1.set_title("Atmospheric Temperature vs. Pressure", fontsize=14)
    ax1.grid(True, which="both", ls=":")
    ax1.legend(fontsize=10)
    
    # Altitude vs Temperature
    ax2.plot(res_planck['T'], res_planck['z'] / 1000, 'r-', lw=2, label="Planck Mean (HITRAN)")
    ax2.plot(res_rosseland['T'], res_rosseland['z'] / 1000, 'b-', lw=2, label="Rosseland Mean (HITRAN)")
    ax2.plot(res_const['T'], res_const['z'] / 1000, 'k--', lw=1.5, label="Constant 1e-4 Baseline")
    
    ax2.axhline(res_planck['z_e']/1000, color='r', linestyle=':', alpha=0.7, label="Planck Emission Height")
    ax2.axhline(res_rosseland['z_e']/1000, color='b', linestyle=':', alpha=0.7, label="Rosseland Emission Height")
    
    ax2.set_xlabel("Temperature (K)", fontsize=12)
    ax2.set_ylabel("Altitude (km)", fontsize=12)
    ax2.set_title("Atmospheric Temperature vs. Altitude", fontsize=14)
    ax2.grid(True, which="both", ls=":")
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    
    plot_path = os.path.join(os.path.dirname(__file__), "temperature_profile.png")
    plt.savefig(plot_path, dpi=300)
    print(f"\nTemperature profile plot saved to '{plot_path}'.")

if __name__ == "__main__":
    run_model()
