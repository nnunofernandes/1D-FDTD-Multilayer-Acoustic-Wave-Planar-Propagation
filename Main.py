# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 17:23:08 2025

@author: nnuno
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd


# This script simulates ultrasound wave propagation using the first-order pressure-velocity formulation
# (momentum and continuity equations) instead of the second-order scalar wave equation (d'Alembert form).

#velocities and rho
c_skin = 1654.93
c_AT = 1461.77
c_cart = 1650.80
c_pzt=4200
c_match=1700
rho_skin = 1131
rho_AT = 927
rho_cart = 1033
rho_pzt=7700
rho_match=1500

freq = 5000e3

# Layer definitions and properties
pzt_thick=c_pzt/(2*freq)
match_thick=c_match/(4*freq)
skin_thick = 0.002
AT_thick = 0.01
cart_thick = 0.0025

# Simulation parameters
L_visible = pzt_thick+match_thick+skin_thick+AT_thick+cart_thick
L_total = 4*L_visible
#Nx = 2400
dx=c_skin/(100*freq)
Nx=int((L_total)/dx)
x = np.linspace(-2*L_visible, 2*L_visible, Nx)
#dx = (L_visible + 0.12) / Nx

# Time step and simulation duration
c_max = c_pzt  # Maximum speed in bone
CFL = 0.5
dt = CFL * dx / c_max

T_period = 1 / freq
T_sim = 0.000012
Nt = int(T_sim / dt)


# Spatial profiles of c and rho
c = np.ones(Nx) * c_cart
rho = np.ones(Nx) * rho_cart

c[(x < pzt_thick+match_thick+skin_thick + AT_thick)] = c_AT
rho[(x < pzt_thick+match_thick+skin_thick + AT_thick)] = rho_AT

c[(x < pzt_thick+match_thick+skin_thick)] = c_skin
rho[(x < pzt_thick+match_thick+skin_thick)] = rho_skin

c[(x < pzt_thick+match_thick)] = c_match
rho[(x < pzt_thick+match_thick)] = rho_match

c[(x < pzt_thick)] = c_pzt
rho[(x < pzt_thick)] = rho_pzt

# Bulk modulus
K = rho * c**2

# Attenuation coefficients (dB/MHz/cm)
alpha_dB_skin = 1.096
alpha_dB_AT = 0.7264
alpha_dB_cart = 4.336
log_to_neper = np.log(10) / 20

alpha_np_skin = alpha_dB_skin * log_to_neper * (freq / 1e6) * 100
alpha_np_AT = alpha_dB_AT * log_to_neper * (freq / 1e6) * 100
alpha_np_cart = alpha_dB_cart * log_to_neper * (freq / 1e6) * 100

# Attenuation profile (in Nepers/m)
alpha = np.ones(Nx) * alpha_np_cart
alpha[(x <pzt_thick+match_thick+skin_thick)] = alpha_np_skin
alpha[(x >= pzt_thick+match_thick+skin_thick) & (x < pzt_thick+match_thick+skin_thick+AT_thick)] = alpha_np_AT
alpha[(x <pzt_thick+match_thick)] = 0

# Damping factor per time step
gamma = alpha * c * dt
gamma = np.clip(gamma, 1e-3, 1)

# Geometric spreading correction (3D)
source_pos = np.argmin(np.abs(x - 0))
r = np.abs(x - x[source_pos]) + dx  # avoid divide by zero
geo_corr = 1 / r**2
geo_corr /= np.max(geo_corr)  # normalize to keep dynamic range reasonable

# Initialize pressure and particle velocity (staggered grid)
p = np.zeros(Nx)
v = np.zeros(Nx)

# Source and receiver positions
skin_transducer_pos=np.argmin(np.abs(x - (pzt_thick+match_thick)))
skin_AT_pos = np.argmin(np.abs(x - (pzt_thick+match_thick+skin_thick)))
AT_cart_pos = np.argmin(np.abs(x - (pzt_thick+match_thick+skin_thick+AT_thick)))
cart_pos = np.argmin(np.abs(x - (pzt_thick+match_thick+skin_thick+AT_thick+0.00125)))
#brain_10_pos = np.argmin(np.abs(x - (pzt_thick+match_thick+skin_thick+bone_thick+0.01)))
#brain_15_pos = np.argmin(np.abs(x - L_visible))
skin_transducer_signal = []
skin_AT_signal = []
AT_cart_signal = []
cart_signal = []
#brain_10_signal = []
#brain_15_signal = []

# Source signal (smooth full sine wave)
source_amplitude = 10
#n_source = int(T_period / dt)
#window = np.hanning(2 * n_source)
#smooth_source = source_amplitude * window[:n_source] * np.sin(2 * np.pi * freq * np.arange(n_source) * dt)

# Plot setup
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 3]})
line, = ax1.plot(x * 1000, np.abs(p) + 1e-12)
ax1.set_xlim(0, L_visible * 1000)
#ax1.set_yscale('log')
ax1.set_ylim(-15, 15)
ax1.set_xlabel("Position [mm]")
ax1.set_ylabel("Pressure Amplitude [Pa]")
ax1.set_title("1D Ultrasound FDTD - Pressure-Velocity Formulation")
time_text = ax1.text(0.02, 0.05, '', transform=ax1.transAxes, fontsize=10, verticalalignment='bottom')

ax1.axvline((pzt_thick+match_thick+skin_thick) * 1000, color='k', linestyle='--', label='Skin-Adipose Tissue Boundary')
ax1.axvline((pzt_thick+match_thick+skin_thick + AT_thick) * 1000, color='r', linestyle='--', label='Adipose Tissue-Cartilage Boundary')
ax1.axvline((pzt_thick+match_thick) * 1000, color='b', linestyle='--', label='Matching Layer-Skin Boundary')
ax1.axvline((pzt_thick) * 1000, color='g', linestyle='--', label='PZT Disc-Matching Layer Boundary')
ax1.legend()

# Time-series plot
time_series_line1, = ax2.plot([], [], label="Transducer-Skin")
time_series_line2, = ax2.plot([], [], label="Skin-Adipose Tissue")
time_series_line3, = ax2.plot([], [], label="Adipose Tissue-Cartilage")
time_series_line4, = ax2.plot([], [], label="Cartilage @ 1.25 mm")
#time_series_line5, = ax2.plot([], [], label="Brain @ 10 mm")
#time_series_line6, = ax2.plot([], [], label="Brain @ 15 mm")
ax2.set_xlim(0, T_sim * 1e6)
#ax2.set_yscale('log')
ax2.set_ylim(-2, 2)
ax2.set_xlabel("Time [µs]")
ax2.set_ylabel("Pressure [Pa]")
ax2.legend()

# Update function
def update(frame):
    global p, v

    # Velocity update (staggered grid)
    v[:-1] -= dt / rho[:-1] * (p[1:] - p[:-1]) / dx

    # Pressure update
    p[1:] -= dt * K[1:] * (v[1:] - v[:-1]) / dx

    # Apply attenuation to pressure field
    p *= (1 - gamma)
    #p *= geo_corr
    
    # Filter low-amplitude noise
    p[np.abs(p) < 1e-7] = 0.0
    #p[1:-1] = 0.25 * p[:-2] + 0.5 * p[1:-1] + 0.25 * p[2:]
    
    # Source injection
    #if frame < len(smooth_source):
        #p[source_pos] = smooth_source[frame]
    if T_period>(frame*dt):
        v[source_pos] = source_amplitude*np.sin(2*np.pi*freq*frame*dt)/rho[source_pos]/c[source_pos]
    elif 2*T_period>(frame*dt):
        v[source_pos] *= 0.1

    # Record data
    skin_transducer_signal.append(p[skin_transducer_pos])
    skin_AT_signal.append(p[skin_AT_pos])
    AT_cart_signal.append(p[AT_cart_pos])
    cart_signal.append(p[cart_pos])
    #brain_10_signal.append(p[brain_10_pos])
    #brain_15_signal.append(p[brain_15_pos])

    # Update plot
    line.set_ydata(p)
    time_text.set_text(f"Time: {frame * dt * 1e6:.2f} µs")


    
    t = np.arange(len(skin_AT_signal)) * dt * 1e6
    time_series_line1.set_data(t, skin_transducer_signal)
    time_series_line2.set_data(t, skin_AT_signal)
    time_series_line3.set_data(t, AT_cart_signal)
    time_series_line4.set_data(t, cart_signal)
    #time_series_line5.set_data(t, brain_10_signal)
    #time_series_line6.set_data(t, brain_15_signal)
    # Save as DataFrame
    if frame==(Nt-1):
        df = pd.DataFrame({
            'Time_us': t,
            'Tranducer-Skin (Pa)': skin_transducer_signal, 
            'Skin-Adipose Tissue (Pa)': skin_AT_signal,
            'Adipose Tissue-Cartilage (Pa)': AT_cart_signal,
            'Cartilage @ 1.25 mm (Pa)': cart_signal})
        df.to_csv('5000k_cartilage.csv', index=False)
    return line, time_text, time_series_line1, time_series_line2,time_series_line3, time_series_line4

ani = animation.FuncAnimation(fig, update, frames=Nt, interval=0.01, blit=True, repeat=False)
plt.tight_layout()
plt.show()
