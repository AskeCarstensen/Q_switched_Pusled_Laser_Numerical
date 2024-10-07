# Equations for Modeling the Passively Q-switched Laser System

This repository contains a numerical model for simulating the behavior of a passively Q-switched diode-pumped solid-state laser system. The purpose of this model is to analyze the evolution of key variables. The equations are solved using a system of coupled differential equations that describe the dynamics of the laser system.

## 1. Population Inversion in the Laser Crystal
The population inversion, ` N_{LA} `, describes the difference in population between the upper and lower laser levels in the laser crystal. The rate equation for $ N_{LA} $ is given by:

$$
\frac{dN_{LA}}{dt} = \frac{P_{pump}}{h \nu_{pump} V_{LA}} - \frac{N_{LA}}{t_{sp}} - 2 N_{LA} \sigma_{LA} \Phi_{LA}
$$

- **Term 1**: $ \frac{P_{pump}}{h \nu_{pump} V_{LA}} $ represents the rate at which pump power $ P_{pump} $ adds population to the upper laser level. Here, $ h $ is Planck's constant, $ \nu_{pump} $ is the pump frequency, and $ V_{LA} $ is the volume of the laser crystal.
- **Term 2**: $ - \frac{N_{LA}}{t_{sp}} $ is the spontaneous emission term, describing how the population decays with the spontaneous emission lifetime $ t_{sp} $.
- **Term 3**: $ - 2 N_{LA} \sigma_{LA} \Phi_{LA} $ is the stimulated emission term, where $ \sigma_{LA} $ is the emission cross-section and $ \Phi_{LA} $ is the circulating photon flux in the laser crystal.

## 2. Circulating Photon Flux in the Laser Crystal
The circulating photon flux in the laser crystal, $ \Phi_{LA} $, describes the number of photons per unit area and time within the laser cavity. The differential equation governing $ \Phi_{LA} $ is:

$$
\frac{d\Phi_{LA}}{dt} = 2 N_{LA} c \sigma_{LA} \frac{\Phi_{LA}}{L_{cav}} - \frac{\Phi_{LA}}{\tau_p} + \frac{N_{LA}}{t_{sp}} \cdot L_{cry} \cdot 10^{-3}
$$

- **Term 1**: $ 2 N_{LA} c \sigma_{LA} \frac{\Phi_{LA}}{L_{cav}} $ describes the amplification of photon flux due to stimulated emission, where $ c $ is the speed of light and $ L_{cav} $ is the roundtrip length of the cavity.
- **Term 2**: $ - \frac{\Phi_{LA}}{\tau_p} $ is the photon loss term, where $ \tau_p $ is the photon lifetime inside the cavity. This lifetime is influenced by the total round-trip losses in the cavity.
- **Term 3**: $ + \frac{N_{LA}}{t_{sp}} \cdot L_{cry} \cdot 10^{-3} $ represents the spontaneous emission contribution to the photon flux.

## 3. Energy Levels of the Saturable Absorber
The saturable absorber is modelled with three main energy levels: $ N_g $ (ground state), $ N_1 $ (intermediate level), and $ N_2 $ (excited state). The rate equations for these levels are:

$$
\frac{dN_g}{dt} = -2 \cdot (N_g - N_2) \sigma_{SA} \Phi_{SA} + \frac{N_1}{t_{1g}}
$$

$$
\frac{dN_1}{dt} = \frac{N_2}{t_{21}} - \frac{N_1}{t_{1g}}
$$

$$
\frac{dN_2}{dt} = 2 \cdot (N_g - N_2) \sigma_{SA} \Phi_{SA} - \frac{N_2}{t_{21}}
$$


## 4. Photon Lifetime in the Cavity
The photon lifetime inside the cavity, $ \tau_p $, is influenced by the total round-trip losses. It is defined as:

$$
\tau_p = \frac{2 L_{cav}}{c} \cdot \frac{1}{L_{roundtrip}}
$$

$$
L_{roundtrip} = L_{passive} + T_{out} + 2 \cdot (N_g - N_2) \sigma_{SA} L_{absorber}
$$


## 5. Output Power
The output power, $ P_{out} $, of the laser is defined as:

$$
P_{out} = T_{out} \cdot \Phi_{LA} \cdot \pi \cdot w_{LA}^2 \cdot \frac{h \cdot c}{\lambda_{LA}}
$$

Where:
- $ T_{out} $ is the transmission of the output mirror.
- $ \Phi_{LA} $ is the circulating photon flux in the laser crystal.
- $ w_{LA} $ is the beam waist in the laser crystal.
- $ \lambda_{LA} $ is the wavelength of the laser emission.

