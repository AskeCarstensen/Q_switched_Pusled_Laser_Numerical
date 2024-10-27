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
sigma_LA = 2.8e-23                      # Laser emission cross section [m^2]ø
sigma_SA = 57e-23                       # Absorber cross section [m^2]
t_sp = 230e-6                           # Spontaneous emission lifetime [s]
t_1g = 3.5e-6                           # Relaxation time [s]
t_21 = 100e-12                          # Relaxation time [s]
c = 3e8                                 # Speed of light [m/s]
h = 6.626e-34                           # Planck's constant [J·s]
wavelength_LA = 1064e-9                 # Laser wavelength [m]
wavelength_pump = 808e-9                # Pump wavelength [m]å
ksi_Cr = 0.0017544                      # Cr concentration in saturable absorber
N0 = 1e26 * ksi_Cr                      # Absorber density [m^-3]
L_pump = 5e-3                           # Length of laser crystal [m]
L_cavity = 200e-3                       # Cavity length [m]
L_absorber = 2e-3                       # Absorber length [m]
T_out = 0.15                            # Output mirror transmission
L_passive = 0.03                        # Passive losses per round trip
w_LA = 150e-6                           # Laser spot size [m]
w_SA = 150e-6                           # New spot size in saturable absorber [m]
A_LA = np.pi * w_LA**2                  # Cross-sectional area of the laser beam [m^2]
E_photon = h * c / wavelength_LA        # Energy per photon [J]
V_LA = L_pump * np.pi * (220e-6)**2     # Volume of laser crystal [m^3]

# Initial conditions
log.info('Setting initial conditions')
N_LA_initial = 0                        # Initial population inversion [m^-3]
Phi_initial = 1e-40                     # Initial photon flux [m^-2·s^-1]
N1_initial = 0                          # Initial ground state population in absorber
N2_initial = 0                          # Initial excited state population in absorber
Ng_initial = N0                         # Initial ground state population in absorber


def laser_system(t, y, P_pump):
    """
    Laser system model equations.

    Parameters:
        - t: Time [s]
        - y: State variables [N_LA, Phi, N_g, N1, N2]
        - P_pump: Pump power [W]
    
    Returns:
        - List of time derivatives for the state variables
    """
    N_LA, Phi, N_g, N1, N2 = y

    # Population inversion equation
    dN_dt = P_pump / (h * c / wavelength_pump * V_LA) - N_LA / t_sp - 2 * N_LA * sigma_LA * Phi

    # Adjust photon flux in saturable absorber based on spot sizes
    Phi_SA = Phi * (w_LA / w_SA)**2

    # Losses from the absober per round trip
    L_SA = 2 * (N_g - N2) * sigma_SA * L_absorber

    # Photon lifetime in the cavity
    t_p = ( 2 * L_cavity / c ) * 1 / (-np.log((1-T_out)*(1-L_passive)*(1-L_SA)))

    # Circulating photon flux equation
    dPhi_dt = 2 * N_LA * c * sigma_LA * Phi * (L_pump/(L_cavity)) - Phi / (t_p) + (N_LA / t_sp) * 1e-3

    # Energy levels in the saturable absorber
    dNg_dt = -2 * (N_g - N2) * sigma_SA * Phi_SA + N1 / t_1g
    dN1_dt = N2 / t_21 - N1 / t_1g
    dN2_dt = 2 * (N_g - N2) * sigma_SA * Phi_SA - N2 / t_21

    return [dN_dt, dPhi_dt, dNg_dt, dN1_dt, dN2_dt]

# Range of pump power
P_pump_values = np.linspace(0.1, 4.0, 10)  
steady_state_results = [] 
pulse_separations = [] 
time_span = (0, 1e-3)

# Solve the system for each pump power and calculate pulse separation
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

    # Extract results for the photon flux
    time_points = np.linspace(time_span[0], time_span[1], 500000)
    results = solution.sol(time_points)
    Phi_results = results[1]

    # Detect pulses based on photon flux
    pulse_times = []  # List to store times of detected pulses
    threshold = np.max(Phi_results) * 0.5  # Define a threshold to detect pulses

    # Detect pulses by checking where photon flux crosses the threshold
    for i in range(1, len(Phi_results) - 1):
        if Phi_results[i - 1] < threshold and Phi_results[i] >= threshold:
            pulse_times.append(time_points[i])

    pulse_intervals = np.diff(pulse_times)  # Time differences between consecutive pulses
    average_pulse_separation = np.mean(pulse_intervals)  # Average time between pulses

    pulse_separations.append(average_pulse_separation)

plt.figure(figsize=(8, 6))
plt.plot(P_pump_values, pulse_separations, marker='o', linestyle='-')
plt.title("Pulse Separation vs Pump Power")
plt.xlabel("Pump Power [W]")
plt.ylabel("Pulse Separation [s]")
plt.grid(True)
plt.show()
