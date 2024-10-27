import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import logging as log

# Assumtions
# we did not consider the electron going from N2 to N1 in the absorber
# We did not consider pump saturation
# We did not consider reffractiv index like self forcusing 
# we did not consider the spatial distribution of the pump beam
# we except a constant field in the cavity
# we did not consider the thermal effects
# we did not consider that the gain medium aborsotion the field
# we did not consider thermal effetcs 

log.basicConfig(
    level=log.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Constants
log.info('Defining constants.')
sigma_LA = 2.8e-23                  # Laser emission cross section [m^2]
sigma_SA = 57e-23                   # Absorber cross section [m^2]
t_sp = 230e-6                       # Spontaneous emission lifetime [s]
t_1g = 3.5e-6                       # Relaxation time [s]
t_21 = 100e-12                      # Relaxation time [s]
c = 3e8                             # Speed of light [m/s]
h = 6.626e-34                       # Planck's constant [J·s]
wavelength_LA = 1064e-9             # Laser wavelength [m]
wavelength_pump = 808e-9            # Pump wavelength [m]
V_LA = 5e-3 * np.pi * (150e-6)**2   # Volume of laser crystal [m^3]
N0 = 1e26                           # Absorber density [m^-3]
L_pump = 5e-3                       # Length of laser crystal [m]
L_cavity = 200e-3                   # Cavity length [m]
L_absorber = 2e-3                   # Absorber length [m]
T_out = 0.15                        # Output mirror transmission
L_passive = 0.03                    # Passive losses per round trip
w_LA = 150e-6                       # Laser spot size [m]a
w_SA = 150e-6                       # New spot size in saturable absorber [m]

# Initial conditions
log.info('Setting initial conditions')
N_LA_initial = 0                    # Initial population inversion [m^-3]
Phi_initial = 1e-40                 # Initial photon flux [m^-2·s^-1]
N1_initial = N0                     # Initial ground state population in absorber
N2_initial = 0                      # Initial excited state population in absorber
Ng_initial = N0                     # Initial ground state population in absorber


def laser_system(t, y, P_pump):
    """
    Laser system model.
    - t: Time [s]
    - y: State variables [N_LA, Phi, N_g, N1, N2]
    - P_pump: Pump power [W]
    
    Returns:
    - List of time derivatives for the state variables
    """
    N_LA, Phi, N_g, N1, N2 = y

    # Population inversion equation
    dN_dt = P_pump / (h * c / wavelength_pump * V_LA) - N_LA / t_sp - 2 * N_LA * sigma_LA * Phi

    # photon lifetime in cavity
    Phi_SA = Phi * (w_SA / w_LA)**2
    L_roundtrip = L_passive + T_out + 2 * (N_g - N2) * sigma_SA  * (L_absorber)
    t_p = (2*L_cavity/c) * 1/L_roundtrip

    # Circulating photon flux equation
    dPhi_dt = 2 * N_LA * c * sigma_LA * Phi / (L_cavity / 2) - Phi / (t_p) + (N_LA / t_sp) * L_pump * 1e-3

    # Energy levels in the saturable absorber
    dNg_dt = -2 * (N_g - N2) * sigma_SA * Phi_SA + N1 / t_1g
    dN1_dt = N2 / t_21 - N1 / t_1g
    dN2_dt = 2 * (N_g - N2) * sigma_SA * Phi_SA - N2 / t_21

    return [dN_dt, dPhi_dt, dNg_dt, dN1_dt, dN2_dt]



# Range of pump power
P_pump_values = np.linspace(4, 4.0, 1)  
steady_state_results = [] 
time_span = (0, 1e-3)

#Photon Flux
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.title("Photon Flux Over Time for Different Pump Powers (Spot Size: {w_SA}} µm)")
plt.xlabel("Time [s]")
plt.ylabel("Photon Flux [$m^{-2}s^{-1}$]")

#Population Inversion
plt.subplot(1, 2, 2)
plt.title("Population Inversion Over Time for Different Pump Powers (Spot Size: {w_SA} µm)")
plt.xlabel("Time [s]")
plt.ylabel("Population Inversion [$m^{-3}$]")

# Sweep over pump powers and solve solve_ivp
for P_pump in P_pump_values:
    log.info(f'Simulating system for pump power: {P_pump:.2f} W')
    initial_conditions = [N_LA_initial, Phi_initial, Ng_initial, N1_initial, N2_initial]
    solution = solve_ivp(
        laser_system, 
        time_span, 
        initial_conditions,
        args=(P_pump,), 
        method='BDF', 
        dense_output=True,
    )

    N_LA_results = solution.y[0]
    Phi_results = solution.y[1]

    # Plot Photon Flux
    plt.subplot(1, 2, 1)
    plt.plot(solution.t, Phi_results, label=f'Pump Power: {P_pump:.2f} W')

    # Plot Population Inversion
    plt.subplot(1, 2, 2)
    plt.plot(solution.t, N_LA_results, label=f'Pump Power: {P_pump:.2f} W')

plt.subplot(1, 2, 1)
plt.legend()
plt.subplot(1, 2, 2)
plt.legend()
plt.tight_layout()
plt.show()