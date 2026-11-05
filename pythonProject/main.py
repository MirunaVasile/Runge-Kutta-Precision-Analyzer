#!/usr/bin/env python3
"""
Simulare numerică RK4 — Pendulul simplu
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 13,
    'legend.fontsize': 10,
    'axes.titlesize': 13,
    'figure.dpi': 150,
})


# PARAMETRI

L     = 1.0
g     = 9.81
omega = np.sqrt(g / L)       # ≈ 3.1321 rad/s
T0    = 2 * np.pi / omega    # ≈ 2.0064 s

print("=" * 60)
print("PARAMETRI SIMULARE")
print(f"  L  = {L} m")
print(f"  g  = {g} m/s²")
print(f"  ω  = {omega:.6f} rad/s")
print(f"  T₀ = {T0:.6f} s  (perioada liniară)")
print("=" * 60)


# ALGORITMUL RK4


def solve_pendulum(theta0_deg, t_end, h=0.005, c=0.0):
    """
    RK4 pentru θ'' + c·θ' + ω²·sin(θ) = 0.
    Condiții inițiale: θ(0)=θ₀, θ'(0)=0.
    Returnează: t, θ(t) [rad], v(t) [rad/s]
    """
    theta0 = np.radians(theta0_deg)

    def F(Y):
        th, v = Y
        return np.array([v, -c * v - omega**2 * np.sin(th)])

    n  = int(t_end / h) + 1
    t  = np.linspace(0, t_end, n)
    Y  = np.zeros((n, 2))
    Y[0] = [theta0, 0.0]

    for i in range(n - 1):
        yn = Y[i]
        k1 = h * F(yn)
        k2 = h * F(yn + k1 / 2)
        k3 = h * F(yn + k2 / 2)
        k4 = h * F(yn + k3)
        Y[i + 1] = yn + (k1 + 2*k2 + 2*k3 + k4) / 6

    return t, Y[:, 0], Y[:, 1]


def compute_period_rk4(theta0_deg, h=0.0005):
    """
    Detectează perioada prin prima traversare a lui v=0
    dinspre negativ spre pozitiv (jumătate de perioadă).
    T = 2 * T_half.
    """
    t_max = max(40.0, 8 * T0)
    t, theta, v = solve_pendulum(theta0_deg, t_max, h=h)

    for i in range(10, len(v) - 1):
        if v[i] < 0 and v[i + 1] >= 0:
            frac   = -v[i] / (v[i + 1] - v[i])
            t_half = t[i] + frac * h
            return 2.0 * t_half

    return None



# FIGURA 1: Comparație liniar vs RK4


fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
t_span = 12.0

cfg = [
    (10,  'steelblue',  r'$\theta_0 = 10°$'),
    (60,  'firebrick',  r'$\theta_0 = 60°$'),
    (90,  'seagreen',   r'$\theta_0 = 90°$'),
]

for ax, (deg, color, lbl) in zip(axes, cfg):
    theta0 = np.radians(deg)

    t_fine    = np.linspace(0, t_span, 5000)
    theta_lin = np.degrees(theta0 * np.cos(omega * t_fine))

    t_rk4, theta_rk4, _ = solve_pendulum(deg, t_span, h=0.005)
    theta_rk4_deg        = np.degrees(theta_rk4)

    T_NL  = compute_period_rk4(deg)
    ratio = T_NL / T0

    ax.plot(t_fine, theta_lin,
            '-',  color=color, lw=2.2,
            label=f'Liniar  —  T = {T0:.4f} s = T₀')
    ax.plot(t_rk4, theta_rk4_deg,
            '--', color=color, lw=2.2, alpha=0.85,
            label=f'RK4  —  T = {T_NL:.4f} s = {ratio:.4f} T₀')

    ax.set_ylabel(r'$\theta(t)$ [grade]')
    ax.set_title(f'{lbl}   (eroare perioadă: {(ratio-1)*100:.2f}%)',
                 fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_ylim(-deg * 1.3, deg * 1.3)

axes[-1].set_xlabel(r'$t$ [s]')
plt.suptitle('Comparație model liniar vs. RK4 — pendul ideal (L = 1 m, g = 9.81 m/s²)',
             fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('fig1_comparatie.pdf', bbox_inches='tight')
plt.savefig('fig1_comparatie.png', bbox_inches='tight', dpi=150)
print("\nSalvat: fig1_comparatie.pdf / .png")
plt.close()


# TABELUL PERIOADELOR

angles_tbl = [5, 10, 20, 30, 45, 60, 90, 120, 150]

print("\n" + "=" * 72)
print(f"  {'θ₀':>5}  |  {'k=sin(θ₀/2)':>12}  |  {'T_RK4 [s]':>10}  |"
      f"  {'T/T₀':>7}  |  {'Eroare':>8}")
print("-" * 72)

table_rows = []
for deg in angles_tbl:
    k     = np.sin(np.radians(deg / 2))
    T_NL  = compute_period_rk4(deg, h=0.0002)
    ratio = T_NL / T0
    err   = (ratio - 1) * 100
    table_rows.append((deg, k, T_NL, ratio, err))
    print(f"  {deg:>4}°  |  {k:>12.6f}  |  {T_NL:>10.6f}  |"
          f"  {ratio:>7.5f}  |  {err:>7.3f}%")

print("=" * 72)
print(f"  T₀ (liniar) = {T0:.6f} s")


# FIGURA 2: T(θ₀)/T₀ vs θ₀

angles_dense = list(range(5, 171, 5))
results      = []
for deg in angles_dense:
    T = compute_period_rk4(deg, h=0.001)
    if T:
        results.append((deg, T / T0))

ang_arr   = np.array([r[0] for r in results])
ratio_arr = np.array([r[1] for r in results])

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(ang_arr, ratio_arr, 'b-o', ms=4, lw=2.2,
        label=r'$T(\theta_0)/T_0$  (RK4)')
ax.axhline(1.0,  color='red',  lw=2.0, ls='--',
           label=r'Aproximare liniară ($T = T_0$)')
ax.axhline(1.02, color='gray', lw=1.2, ls=':',
           label='Toleranță 2%')

ax.fill_between([0, 20], [1.0, 1.0], [1.02, 1.02],
                color='green', alpha=0.15)
ax.text(10, 1.005, 'valid', color='darkgreen',
        ha='center', fontsize=10, fontweight='bold')

# Annotate key points from table
for deg, k, T_NL, ratio, err in table_rows:
    if deg in [30, 60, 90, 120]:
        ax.plot(deg, ratio, 'bo', ms=8)
        ax.annotate(f'{ratio:.3f}',
                    (deg, ratio),
                    textcoords='offset points',
                    xytext=(6, 5), fontsize=9, color='navy')

ax.set_xlabel(r'Amplitudinea inițială $\theta_0$ [grade]')
ax.set_ylabel(r'$T(\theta_0) \;/\; T_0$')
ax.set_xlim(0, 172)
ax.set_ylim(0.98, ratio_arr.max() * 1.06)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_title(r'Raportul de perioade $T(\theta_0)/T_0$ — calculat prin RK4')
plt.tight_layout()
plt.savefig('fig2_perioada.pdf', bbox_inches='tight')
plt.savefig('fig2_perioada.png', bbox_inches='tight', dpi=150)
print("\nSalvat: fig2_perioada.pdf / .png")
plt.close()


# FIGURA 3: Disiparea energiei (pendul realist, RK4)

theta0_deg = 30.0
t_end      = 15.0

fig, ax = plt.subplots(figsize=(11, 6))

damp_cfg = [
    (0.2, 'steelblue', r'$c = 0.2\;\mathrm{s}^{-1}$  ($\tau = 5$ s)'),
    (0.5, 'firebrick', r'$c = 0.5\;\mathrm{s}^{-1}$  ($\tau = 2$ s)'),
    (1.0, 'seagreen',  r'$c = 1.0\;\mathrm{s}^{-1}$  ($\tau = 1$ s)'),
]

for c_val, color, label in damp_cfg:
    t, theta, v = solve_pendulum(theta0_deg, t_end, h=0.005, c=c_val)
    E  = v**2 / 2.0 + omega**2 * (1.0 - np.cos(theta))
    E0 = E[0]
    ax.plot(t, E / E0, '-', color=color, lw=2.5, label=f'RK4: {label}')

    tau = 1.0 / c_val
    ax.axvline(tau, color=color, lw=1.0, ls=':', alpha=0.7)

t_fine = np.linspace(0, t_end, 2000)
ax.plot(t_fine, np.exp(-0.5 * t_fine),
        ':', color='firebrick', lw=1.8, alpha=0.7,
        label=r'Anvelopă $e^{-ct}$, $c=0.5$ (predicție liniară)')

ax.axhline(np.exp(-1), color='gray', lw=1.5, ls='--')
ax.text(13.0, np.exp(-1) + 0.025,
        r'$e^{-1} \approx 0.368$', color='gray', fontsize=10)

ax.set_xlabel(r'$t$ [s]')
ax.set_ylabel(r'$E(t) \;/\; E(0)$')
ax.set_xlim(0, t_end)
ax.set_ylim(-0.02, 1.08)
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_title(r'Disiparea energiei mecanice — pendul realist RK4 ($\theta_0 = 30°$)')
plt.tight_layout()
plt.savefig('fig3_disipare.pdf', bbox_inches='tight')
plt.savefig('fig3_disipare.png', bbox_inches='tight', dpi=150)
print("Salvat: fig3_disipare.pdf / .png")
plt.close()


# TABEL (valori din RK4)

print("\n\n=== TABEL LaTeX ===\n")
print(r"\begin{table}[H]")
print(r"\centering")
print(r"\caption{Perioada pendulului neliniar calculat\u{a} numeric prin RK4,")
print(r"         pentru $L=1\,\mathrm{m}$, $g=9{,}81\,\mathrm{m/s}^2$,")
print(f"         $T_0 = {T0:.4f}\\,\\mathrm{{s}}$, pas $h=0{{,}}0002\\,\\mathrm{{s}}$.}}")
print(r"\label{tab:perioade_rk4}")
print(r"\renewcommand{\arraystretch}{1.3}")
print(r"\begin{tabular}{ccccc}")
print(r"\hline")
print(r"$\theta_0$ & $k = \sin(\theta_0/2)$ & "
      r"$T_{\mathrm{RK4}}/T_0$ & $T_{\mathrm{RK4}}$ [s] & Eroare relativ\u{a} \\")
print(r"\hline")
for deg, k, T_NL, ratio, err in table_rows:
    k_str   = f"{k:.4f}".replace('.', '{,}')
    r_str   = f"{ratio:.4f}".replace('.', '{,}')
    T_str   = f"{T_NL:.4f}".replace('.', '{,}')
    e_str   = f"{err:.3f}".replace('.', '{,}')
    print(f"${deg}^\\circ$ & ${k_str}$ & ${r_str}$ & ${T_str}$ & ${e_str}\\%$ \\\\")
print(r"\hline")
print(r"\end{tabular}")
print(r"\end{table}")
print()