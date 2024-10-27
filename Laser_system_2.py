import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import logging as log

log.basicConfig(
    level=log.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Constants
log.info('Defining constants.')
sigma_LA = 2.8e-23                      # Laser emission cross section [m^2]
sigma_SA = 57e-23                       # Absorber cross section [m^2]
t_sp = 230e-6                           # Spontaneous emission lifetime [s]
t_1g = 3.5e-6                           # Relaxation time [s]
t_21 = 100e-12                          # Relaxation time [s]
c = 3e8                                 # Speed of light [m/s]
h = 6.626e-34                           # Planck's constant [J·s]
wavelength_LA = 1064e-9                 # Laser wavelength [m]
wavelength_pump = 808e-9                # Pump wavelength [m]
ksi_Cr = 0.0017544                      # Cr concentration in saturable absorber
N0 = 1e26 * ksi_Cr                      # Absorber density [m^-3]
L_pump = 5e-3                           # Length of laser crystal [m]
L_cavity = 200e-3                       # Cavity length [m]
L_absorber = 2e-3                       # Absorber length [m]
T_out = 0.15                            # Output mirror transmission
L_passive = 0.03                        # Passive losses per round trip
w_LA = 150e-6                           # Laser spot size [m]
V_LA = L_pump * np.pi * (220e-6)**2     # Volume of laser crystal [m^3]

# Initial conditions
log.info('Setting initial conditions')
N_LA_initial = 0                        # Initial population inversion [m^-3]
Phi_initial = 1e-40                     # Initial photon flux [m^-2·s^-1]
N1_initial = 0                          # Initial ground state population in absorber
N2_initial = 0                          # Initial excited state population in absorber
Ng_initial = N0                         # Initial ground state population in absorber

# Laser system model equations
def laser_system(t, y, P_pump, w_SA):
    N_LA, Phi, N_g, N1, N2 = y

    dN_dt = P_pump / (h * c / wavelength_pump * V_LA) - N_LA / t_sp - 2 * N_LA * sigma_LA * Phi
    Phi_SA = Phi * (w_LA / w_SA)**2
    L_SA = 2 * (N_g - N2) * sigma_SA * L_absorber
    t_p = (2 * L_cavity / c) * 1 / (-np.log((1 - T_out) * (1 - L_passive) * (1 - L_SA)))
    dPhi_dt = 2 * N_LA * c * sigma_LA * Phi * (L_pump / L_cavity) - Phi / t_p + (N_LA / t_sp) * 1e-3
    dNg_dt = -2 * (N_g - N2) * sigma_SA * Phi_SA + N1 / t_1g
    dN1_dt = N2 / t_21 - N1 / t_1g
    dN2_dt = 2 * (N_g - N2) * sigma_SA * Phi_SA - N2 / t_21

    return [dN_dt, dPhi_dt, dNg_dt, dN1_dt, dN2_dt]

# Settings
P_pump = 4.0  # Fixed pump power of 4 W
w_SA_values = [100e-6, 150e-6]  # Spot sizes to test
time_span = (0, 1e-3)

# Plotting setup
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.title("Photon Flux Over Time for w_SA = 100 µm and 150 µm")
plt.xlabel("Time [s]")
plt.ylabel("Photon Flux [$m^{-2}s^{-1}$]")

plt.subplot(1, 2, 2)
plt.title("Population Inversion Over Time for w_SA = 100 µm and 150 µm")
plt.xlabel("Time [s]")
plt.ylabel("Population Inversion [$m^{-3}$]")

# Solve for each w_SA value
for w_SA in w_SA_values:
    log.info(f'Simulating system for spot size w_SA: {w_SA*1e6:.0f} µm')
    initial_conditions = [N_LA_initial, Phi_initial, Ng_initial, N1_initial, N2_initial]
    solution = solve_ivp(
        laser_system, 
        time_span, 
        initial_conditions,
        args=(P_pump, w_SA), 
        method='BDF', 
        dense_output=True,
    )
    
    N_LA_results = solution.y[0]
    Phi_results = solution.y[1]

    # Plot Photon Flux
    plt.subplot(1, 2, 1)
    plt.plot(solution.t, Phi_results, label=f'w_SA: {w_SA*1e6:.0f} µm')

    # Plot Population Inversion
    plt.subplot(1, 2, 2)
    plt.plot(solution.t, N_LA_results, label=f'w_SA: {w_SA*1e6:.0f} µm')

plt.subplot(1, 2, 1)
plt.legend()
plt.subplot(1, 2, 2)
plt.legend()
plt.tight_layout()
plt.show()
