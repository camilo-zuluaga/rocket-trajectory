#####################################################
## Simple Rocket Simulation using ODEs into Orbit  ##
#####################################################


"""
    Imports
    Matplotlib: Plotting library
    Numpy: Numerical Python library
    Scipy: Scientific Python library for solving differential equations
"""

import matplotlib.pyplot as plt
import scienceplots
import numpy as np
from scipy.integrate import solve_ivp

plt.style.use("dark_background")
plt.style.use(["science", "no-latex"])

# ------------------------- Constants -------------------------
G0 = 9.81
DEG = np.pi / 180
RHO0 = 1.225
Re = 6371000
HSCALE = 7500
OMEGA = 7.292115e-5
CD = 0.74

# ------------------------ Rocket inputs ------------------------

# Rocket dimensions
diam = 1.2
A = np.pi / 4 * (diam) ** 2

# Rocket structure and propellant mass
mprop = 10000  # KG
mprop2 = 2000  # KG

mstruc = 900  # KG
mstruc2 = 300  # KG

# Payload mass
mpl = 300  # KG

# Liftoff mass
m0 = mprop + mprop2 + mstruc + mstruc2 + mpl
m0s2 = mprop2 + mstruc2 + mpl

# Thrust Burn times
tburn = 138.69
tburn2 = 252.89


# Flow rates
m_dot = mprop / tburn
m_dot2 = mprop2 / tburn2

# Thrust values
Thrust = 218.985 * 1000
Thrust2 = 28.192 * 1000

# Pitch over altitude
hturn = 1430

# Time range for simulation
max_tburn = tburn + tburn2


# ------------------------- Initial conditions for ODES -------------------------
t_max = max_tburn
v0 = 0
psi0 = 10.489 * DEG
theta0 = 0
h0 = 0

# ------------------------- Arrays for post processing -------------------------
mass = []
dynamic_pressure = []
thrust_force = []
drag_force = []
gravity_force = []


# ------------------------- Differential equations -------------------------
def derivatives(t, y):
    v = y[0]
    psi = y[1]
    theta = y[2]
    h = y[3]

    g = G0 / (1 + h / Re) ** 2
    rho = RHO0 * np.exp(-h / HSCALE)
    D = 1 / 2 * rho * v**2 * A * CD
    q = 1 / 2 * rho * v**2 / 1000

    drag_force.append(D)
    dynamic_pressure.append(q)

    if t < tburn:
        dm_dt = m0 - m_dot * t
        T = Thrust
    elif t < tburn + tburn2:
        dm_dt = m0s2 - m_dot2 * (t - tburn)
        T = Thrust2
    else:
        dm_dt = m0s2 - mprop2
        T = 0

    gravity_force.append(dm_dt * g)
    thrust_force.append(T)
    mass.append(dm_dt)

    if h <= hturn:
        psi_dot = 0
        v_dot = T / dm_dt - D / dm_dt - g
        theta_dot = 0
        h_dot = v
    else:
        phi_dot = g * np.sin(psi) / v
        v_dot = T / dm_dt - D / dm_dt - g * np.cos(psi)
        h_dot = v * np.cos(psi)
        theta_dot = v * np.sin(psi) / (Re + h)
        psi_dot = phi_dot - theta_dot
    return [v_dot, psi_dot, theta_dot, h_dot]


# ------------------------- Solve ODE System --------------------------------
sol = solve_ivp(derivatives, [0, t_max], [v0, psi0, theta0, h0], max_step=1)

# ------------------------- Post processing results -------------------------
vrel = sol.y[0] / 1000
vabs = vrel + Re * OMEGA / 1000
psi = sol.y[1]
psideg = psi / DEG
theta = sol.y[2]
dr = theta * Re / 1000
h = sol.y[3] / 1000
htot = h + Re / 1000
rearraytheta = np.linspace(0, 2 * np.pi, 100)
rearray = np.full((100, 1), Re / 1000)
tS = sol.t
tA = np.linspace(0, max_tburn, len(thrust_force))
tD = np.linspace(0, h[-1], len(thrust_force))


# ------------------------- Plotting -------------------------

# ~~~~~~~~~~~~~~~~~~~~~~ Plotting height vs downrange ~~~~~~~~~~~~~~~~~~~~~~
# plt.plot(dr, h, zorder=3, color="#0080FF")
# plt.ylabel("$Height\ (km)$")
# plt.xlabel("$Downrange\ (km)$")
# plt.axhline(y=100, color="#CBCBCB", zorder=2)
# plt.text(
#     x=1.03 * max(dr),
#     y=100 + 5,
#     s="$Karman\ Line$",
#     horizontalalignment="right",
#     verticalalignment="bottom",
#     color="#A4A4A4",
#     zorder=4,
# )
# stageSep = plt.scatter(102.45, 67.8404088, color="#FF4600", zorder=4, alpha=0.9)
# orbitReach = plt.scatter(540, 160, color="#7C00FF", zorder=4, alpha=0.9)
# altitude = plt.scatter(1200, 184, color="#351FEC", zorder=4, alpha=0.9)
# plt.legend(
#     [stageSep, orbitReach, altitude],
#     ["Stage Separation", "Orbit Reach", "184 Km"],
#     fontsize=8,
#     loc="upper left",
# )
# plt.axhline(y=0, color="#008C13")
# plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1)
# plt.savefig("height_vs_downrange.png", dpi=400)
# plt.show()

# # ~~~~~~~~~~~~~~~~~~~~~~ Plotting Forces vs Time ~~~~~~~~~~~~~~~~~~~~~~

# plt.plot(tA, np.array(thrust_force) / 1000, color="#0080FF")
# plt.plot(tA, np.array(drag_force) / 1000, color="#ecb120")
# plt.plot(tA, np.array(gravity_force) / 1000, color="#d8521a")
# plt.ylabel("$Forces\ (Kilo Newtons)$")
# plt.xlabel("$Time\ (s)$")
# plt.legend(["Thrust", "Drag", "Gravity"], loc="upper right")
# plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1)
# plt.savefig("forces_vs_time.png", dpi=400)
# plt.show()

# # ~~~~~~~~~~~~~~~~~~~~~~ Plotting Mass vs Time ~~~~~~~~~~~~~~~~~~~~~~

# plt.plot(tA, mass, color="#0080FF")
# plt.ylabel("$Mass\ (Kg)$")
# plt.xlabel("$Time\ (s)$")
# stage_sep = plt.scatter(tburn, 3650, color="#FF4600", zorder=4, alpha=0.9)
# plt.legend(
#     [stage_sep],
#     ["Stage Separation"],
#     fontsize=8,
#     loc="upper right",
# )
# plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1)
# plt.savefig("mass_vs_time.png", dpi=400)
# plt.show()

# # ~~~~~~~~~~~~~~~~~~~~~~ Plotting Velocity vs Time ~~~~~~~~~~~~~~~~~~~~~~

# plt.plot(tS, vrel, color="#0080FF")
# plt.ylabel("$Velocity\ (Km/s)$")
# plt.xlabel("$Time\ (s)$")
# stage_sep = plt.scatter(tburn, 2.93, color="#FF4600", zorder=4, alpha=0.9)
# orbit_reached = plt.scatter(390.247, 7.8, color="#7e2e8c", zorder=4, alpha=0.9)
# plt.legend(
#     [stage_sep, orbit_reached], ["Stage Separation", "Orbit Speed"], loc="upper left"
# )
# plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1)
# plt.savefig("velocity_vs_time.png", dpi=400)
# plt.show()

# # ~~~~~~~~~~~~~~~~~~~~~~ Plotting Dynamic Pressure vs Height ~~~~~~~~~~~~~~~~~~~~~~

# plt.plot(tD, dynamic_pressure, color="#0080FF")
# plt.ylabel("$Dynamic Pressure\ (KPa)$")
# plt.xlabel("$Height\ (Km)$")
# plt.text(
#     x=190,
#     y=30.5,
#     s="$Maximum\ Dynamic\ Pressure$",
#     horizontalalignment="right",
#     verticalalignment="bottom",
#     color="#d8521a",
# )
# plt.axhline(y=33.5, color="#d8521a", linestyle="--")
# plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1)
# plt.savefig("dynamic_pressure_vs_height.png", dpi=400)
# plt.show()

# # ~~~~~~~~~~~~~~~~~~~~~~ Plotting Altitude vs Velocity ~~~~~~~~~~~~~~~~~~~~~~

# plt.plot(vrel, h, color="#0080FF")
# plt.ylabel("$Height\ (Km)$")
# plt.xlabel("$Velocity\ (Km/s)$")
# stage_sep = plt.scatter(2.95, 68.37, color="#FF4600", zorder=4, alpha=0.9)
# orbit_reached = plt.scatter(7.85, 184.2, color="#7e2e8c", zorder=4, alpha=0.9)
# plt.legend(
#     [stage_sep, orbit_reached], ["Stage Separation", "Orbital Speed"], loc="upper left"
# )
# plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1)
# plt.savefig("altitude_vs_velocity.png", dpi=400)
# plt.show()


# ~~~~~~~~~~~~~~~~~~~~~~ Plotting all in one ~~~~~~~~~~~~~~~~~~~~~~
# Comment all of this and uncomment from line 149 to 234 to see individual plots
plt.figure(1)
# Left side plots
plt.subplot(231)
plt.plot(dr, h, zorder=3, color="#0080FF")
plt.ylabel("$Height\ (km)$")
plt.xlabel("$Downrange\ (km)$")
plt.axhline(y=100, color="#CBCBCB", zorder=2)
plt.text(
    x=1.03 * max(dr),
    y=100 + 5,
    s="$Karman\ Line$",
    horizontalalignment="right",
    verticalalignment="bottom",
    color="#A4A4A4",
    zorder=4,
)
stageSep = plt.scatter(109.25, 67.8404088, color="#FF4600", zorder=4, alpha=0.7)
orbitReach = plt.scatter(540, 160, color="#7C00FF", zorder=4, alpha=0.7)
altitude = plt.scatter(1200, 184, color="#351FEC", zorder=4, alpha=0.7)
plt.legend(
    [stageSep, orbitReach, altitude],
    ["Stage Separation", "Orbit Reach", "184 Km"],
    fontsize=8,
    loc="upper left",
)
plt.axhline(y=0, color="#008C13")
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)

plt.subplot(232)
plt.plot(tA, np.array(thrust_force) / 1000, color="#0080FF")
plt.plot(tA, np.array(drag_force) / 1000, color="#ecb120")
plt.plot(tA, np.array(gravity_force) / 1000, color="#d8521a")
plt.ylabel("$Forces\ (Kilo Newtons)$")
plt.xlabel("$Time\ (s)$")
plt.legend(["Thrust", "Drag", "Gravity"], loc="upper right")
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)

plt.subplot(233)
plt.plot(tA, mass, color="#0080FF")
plt.ylabel("$Mass\ (Kg)$")
plt.xlabel("$Time\ (s)$")
stage_sep = plt.scatter(tburn, 3650, color="#FF4600", zorder=4, alpha=0.7)
plt.legend(
    [stage_sep],
    ["Stage Separation"],
    fontsize=8,
    loc="upper right",
)
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)

# Right side plots
plt.subplot(234)
plt.plot(tS, vrel, color="#0080FF")
plt.ylabel("$Velocity\ (Km/s)$")
plt.xlabel("$Time\ (s)$")
stage_sep = plt.scatter(tburn, 2.93, color="#FF4600", zorder=4, alpha=0.7)
orbit_reached = plt.scatter(390.247, 7.8, color="#7e2e8c", zorder=4, alpha=0.7)
plt.legend(
    [stage_sep, orbit_reached], ["Stage Separation", "Orbit Speed"], loc="upper left"
)
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)

plt.subplot(235)
plt.plot(tD, dynamic_pressure, color="#0080FF")
plt.ylabel("$Dynamic Pressure\ (KPa)$")
plt.xlabel("$Height\ (Km)$")
plt.text(
    x=190,
    y=30.5,
    s="$Maximum\ Dynamic\ Pressure$",
    horizontalalignment="right",
    verticalalignment="bottom",
    color="#d8521a",
)
plt.axhline(y=33.5, color="#d8521a", linestyle="--")
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)

plt.subplot(236)
plt.plot(vrel, h, color="#0080FF")
plt.ylabel("$Height\ (Km)$")
plt.xlabel("$Velocity\ (Km/s)$")
stage_sep = plt.scatter(2.95, 68.37, color="#FF4600", zorder=4, alpha=0.9)
orbit_reached = plt.scatter(7.85, 184.2, color="#7e2e8c", zorder=4, alpha=0.9)
plt.legend(
    [stage_sep, orbit_reached], ["Stage Separation", "Orbital Speed"], loc="upper left"
)
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)

plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~ Plotting Trajectory on Earth ~~~~~~~~~~~~~~~~~~~~~~

fig, ax = plt.subplots()

circle = plt.Circle(
    (0, 0),
    6371,
    edgecolor="#0032FF",
    fill=True,
    facecolor="#0071E1",
    alpha=0.4,
    zorder=2,
)
ax.add_artist(circle)
x_traj1 = rearray * np.cos(rearraytheta)
y_traj1 = rearray * np.sin(rearraytheta)
x_traj2 = htot * np.cos(theta)
y_traj2 = htot * np.sin(theta)
ax.plot(y_traj1, x_traj1, color="#5B78EC", lw=1)
ax.plot(y_traj2, x_traj2, color="#5B78EC", lw=1)
# ax.set_xlim([-1000, 2500])
# ax.set_ylim([5500, 6800])
ax.set_aspect("equal")
plt.grid(True, linestyle="--", linewidth=1, color="#ECECEC", zorder=1, alpha=0.2)
plt.show()
