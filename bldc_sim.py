"""
BLDC Motor — Six-Step Trapezoidal Control Simulation with PI Speed Loop
=========================================================================
Models a Y-connected BLDC motor driven by ideal six-step (120-degree
conduction) commutation, with a PI outer speed loop. This is the classic
first-principles control scheme used before moving to sinusoidal FOC.

Model simplifications (documented, not hidden — good talking points in
an interview):
  - Trapezoidal back-EMF approximated via a clipped triangle wave.
  - Two active phases modeled as a single "line" RL circuit (standard
    six-step simplification: i_positive = -i_negative, floating phase = 0).
  - No inner current loop — duty cycle commanded directly by the speed
    PI. A real drive would cascade an inner current loop under this.
  - No PWM switching ripple modeled at the power-electronics level;
    "duty" acts as an averaged applied voltage fraction.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------
# Motor & drive parameters (small BLDC, e.g. drone/robotics scale)
# ---------------------------------------------------------------------
R      = 0.5        # phase resistance [ohm]
L      = 0.5e-3      # phase inductance [H]
Ke     = 0.01         # back-EMF constant [V*s/rad, electrical]
P      = 4           # pole pairs
J      = 1.0e-4       # rotor inertia [kg*m^2]
B_fric = 1.0e-4       # viscous friction [N*m*s]
Vdc    = 24.0         # DC bus voltage [V]

# Sanity check: back-EMF must stay below Vdc across the target speed range,
# or the drive saturates and can never reach the reference. Max mechanical
# speed achievable at duty=1 is roughly Vdc / (Ke * P).
_max_speed_check = Vdc / (Ke * P)

R_line = 2 * R        # two phases in series during conduction
L_line = 2 * L

# PI speed controller gains (hand-tuned below)
Kp = 0.06
Ki = 1.2

# ---------------------------------------------------------------------
# Reference speed profile: step to 300 rad/s, then a load torque step
# ---------------------------------------------------------------------
OMEGA_REF = 300.0      # rad/s (mechanical)
T_LOAD_STEP = 0.15     # s
T_LOAD_NM = 0.02        # N*m

def ref_speed(t):
    return OMEGA_REF

def load_torque(t):
    return T_LOAD_NM if t > T_LOAD_STEP else 0.0

# ---------------------------------------------------------------------
# Trapezoidal back-EMF shape: clipped triangle wave, +-1 amplitude
# ---------------------------------------------------------------------
def trap_shape(theta):
    theta = np.mod(theta, 2 * np.pi)
    tri = 2.0 * (2.0 * np.abs(theta / (2 * np.pi) - np.floor(theta / (2 * np.pi) + 0.5)) - 0.5)
    return np.clip(1.5 * tri, -1.0, 1.0)

def back_emfs(theta_e, omega_mech):
    omega_e = P * omega_mech
    ea = Ke * omega_e * trap_shape(theta_e)
    eb = Ke * omega_e * trap_shape(theta_e - 2 * np.pi / 3)
    ec = Ke * omega_e * trap_shape(theta_e - 4 * np.pi / 3)
    return ea, eb, ec

# ---------------------------------------------------------------------
# Six-step commutation table: which phase is +, -, floating, per sector
# Sector k covers electrical angle [k*60deg, (k+1)*60deg)
# ---------------------------------------------------------------------
# entries: (positive_phase, negative_phase, floating_phase) as indices 0=A,1=B,2=C
SECTOR_TABLE = [
    (0, 1, 2),  # 0-60:   A+ B-
    (0, 2, 1),  # 60-120: A+ C-
    (1, 2, 0),  # 120-180: B+ C-
    (1, 0, 2),  # 180-240: B+ A-
    (2, 0, 1),  # 240-300: C+ A-
    (2, 1, 0),  # 300-360: C+ B-
]

# Commutation must be aligned to where each phase's back-EMF is actually
# flat, not to an arbitrary theta_e=0 reference. trap_shape(theta) for
# phase A is flat +1 over [150deg,210deg) and flat -1 over [330deg,30deg) -
# so the sector table (built assuming sector0 = A's positive flat top)
# needs a 150 degree offset applied before the lookup.
SECTOR_OFFSET = np.deg2rad(150)

def get_sector(theta_e):
    theta_shifted = np.mod(theta_e - SECTOR_OFFSET, 2 * np.pi)
    return int(theta_shifted // (np.pi / 3)) % 6

# ---------------------------------------------------------------------
# State vector: [theta_mech, omega_mech, i_line, xi (speed error integral)]
# ---------------------------------------------------------------------
def deriv(t, y):
    theta_m, omega, i, xi = y
    theta_e = P * theta_m
    sector = get_sector(theta_e)
    pos, neg, flt = SECTOR_TABLE[sector]

    ea, eb, ec = back_emfs(theta_e, omega)
    e_phases = [ea, eb, ec]
    e_diff = e_phases[pos] - e_phases[neg]

    # --- PI speed controller (with anti-windup) ---
    error = ref_speed(t) - omega
    duty_unclamped = Kp * error + Ki * xi
    duty = np.clip(duty_unclamped, -1.0, 1.0)
    v_applied = duty * Vdc
    # anti-windup: stop integrating further in the saturating direction
    saturated = (duty_unclamped != duty)
    dxi_dt = 0.0 if (saturated and np.sign(error) == np.sign(duty_unclamped)) else error

    # --- electrical dynamics (line RL circuit during conduction) ---
    di_dt = (v_applied - i * R_line - e_diff) / L_line

    # --- torque from instantaneous power balance ---
    ia = i if pos == 0 else (-i if neg == 0 else 0.0)
    ib = i if pos == 1 else (-i if neg == 1 else 0.0)
    ic = i if pos == 2 else (-i if neg == 2 else 0.0)
    if omega > 1.0:
        Te = (ea * ia + eb * ib + ec * ic) / omega
    else:
        Te = P * Ke * i  # low-speed approximation to avoid divide-by-near-zero

    domega_dt = (Te - B_fric * omega - load_torque(t)) / J
    dtheta_dt = omega

    return [dtheta_dt, domega_dt, di_dt, dxi_dt]

# ---------------------------------------------------------------------
# Simulate
# ---------------------------------------------------------------------
t_span = (0, 0.3)
t_eval = np.linspace(*t_span, 6000)
y0 = [0.0, 0.0, 0.0, 0.0]

sol = solve_ivp(deriv, t_span, y0, t_eval=t_eval, method="RK45", max_step=2e-5)

theta_m, omega, i_line, xi = sol.y
t = sol.t

# reconstruct per-phase currents & torque for plotting
ia_arr, ib_arr, ic_arr, Te_arr = [], [], [], []
for k in range(len(t)):
    theta_e = P * theta_m[k]
    sector = get_sector(theta_e)
    pos, neg, flt = SECTOR_TABLE[sector]
    ea, eb, ec = back_emfs(theta_e, omega[k])
    e_phases = [ea, eb, ec]
    i = i_line[k]
    ia = i if pos == 0 else (-i if neg == 0 else 0.0)
    ib = i if pos == 1 else (-i if neg == 1 else 0.0)
    ic = i if pos == 2 else (-i if neg == 2 else 0.0)
    ia_arr.append(ia); ib_arr.append(ib); ic_arr.append(ic)
    if omega[k] > 1.0:
        Te = (ea * ia + eb * ib + ec * ic) / omega[k]
    else:
        Te = P * Ke * i
    Te_arr.append(Te)

ia_arr = np.array(ia_arr); ib_arr = np.array(ib_arr); ic_arr = np.array(ic_arr)
Te_arr = np.array(Te_arr)

# ---------------------------------------------------------------------
# Performance metrics (for resume / report)
# ---------------------------------------------------------------------
final_speed = omega[-1]
target_90 = 0.9 * OMEGA_REF
rise_idx = np.argmax(omega >= target_90)
rise_time = t[rise_idx] if omega[rise_idx] >= target_90 else None

overshoot_pct = (np.max(omega) - OMEGA_REF) / OMEGA_REF * 100

# settling time: within 2% band, before load step
band = 0.02 * OMEGA_REF
pre_load_mask = t < T_LOAD_STEP
settle_idx = None
for k in range(len(t)):
    if t[k] > 0.02 and np.all(np.abs(omega[k:np.searchsorted(t, T_LOAD_STEP)] - OMEGA_REF) < band):
        settle_idx = k
        break
settling_time = t[settle_idx] if settle_idx is not None else None

# load-step recovery
load_idx = np.searchsorted(t, T_LOAD_STEP)
dip = OMEGA_REF - np.min(omega[load_idx:load_idx + 500])
recovery_idx = None
for k in range(load_idx, len(t)):
    if np.abs(omega[k] - OMEGA_REF) < band:
        recovery_idx = k
        break
recovery_time = (t[recovery_idx] - T_LOAD_STEP) if recovery_idx is not None else None

print("=== Performance metrics ===")
print(f"Rise time to 90% of {OMEGA_REF} rad/s: {rise_time*1000:.2f} ms" if rise_time else "Rise time: not reached")
print(f"Overshoot: {overshoot_pct:.2f} %")
print(f"Settling time (2% band): {settling_time*1000:.2f} ms" if settling_time else "Settling time: not reached")
print(f"Speed dip on {T_LOAD_NM} N*m load step: {dip:.2f} rad/s")
print(f"Recovery time after load step: {recovery_time*1000:.2f} ms" if recovery_time else "Recovery: not reached")
print(f"Final steady-state speed: {final_speed:.2f} rad/s (target {OMEGA_REF})")

# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(9, 9))

axs[0].plot(t, omega, label="Actual speed", color="#1f77b4")
axs[0].axhline(OMEGA_REF, color="gray", linestyle="--", label="Reference")
axs[0].axvline(T_LOAD_STEP, color="red", linestyle=":", label="Load step")
axs[0].set_ylabel("Speed [rad/s]")
axs[0].set_title("BLDC Six-Step Control: Speed Response")
axs[0].legend(loc="lower right")
axs[0].grid(alpha=0.3)
axs[0].set_xlim(t_span[0], t_span[1])

axs[1].plot(t, ia_arr, label="Phase A", linewidth=0.8)
axs[1].plot(t, ib_arr, label="Phase B", linewidth=0.8, alpha=0.7)
axs[1].plot(t, ic_arr, label="Phase C", linewidth=0.8, alpha=0.7)
axs[1].set_ylabel("Current [A]")
axs[1].set_title("Trapezoidal Phase Currents (Six-Step Commutation)")
axs[1].legend(loc="upper right")
axs[1].grid(alpha=0.3)
axs[1].set_xlim(0.05, 0.08)  # zoom in to actually see commutation steps

axs[2].plot(t, Te_arr, color="#d62728", linewidth=0.8)
axs[2].axvline(T_LOAD_STEP, color="red", linestyle=":")
axs[2].set_ylabel("Torque [N*m]")
axs[2].set_xlabel("Time [s]")
axs[2].set_title("Electromagnetic Torque (note commutation ripple)")
axs[2].grid(alpha=0.3)
axs[2].set_xlim(t_span[0], t_span[1])

plt.tight_layout()
plt.savefig("bldc_sim_results.png", dpi=150)
print("\nSaved plot to bldc_sim_results.png")
