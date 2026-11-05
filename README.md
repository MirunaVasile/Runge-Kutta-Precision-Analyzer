# RK4 Numerical Simulation of the Simple Pendulum

A numerical simulation of the simple pendulum (ideal and damped) written in Python using the **4th-order Runge-Kutta method (RK4)**. The project compares the **linear** model (`sin θ ≈ θ`) with the full **nonlinear** model and determines the range in which the linear approximation remains valid.

The code is part of a bachelor's thesis in mathematics, **Chapter 4 – Mathematical Models Governed by Differential Equations** (section 4.3, *Numerical simulation: comparison between the linear and nonlinear models*).

---

## What the project does

- Integrates the pendulum equation with RK4, implemented from scratch (no off-the-shelf ODE solvers).
- Compares the linear solution `θ(t) = θ₀·cos(ωt)` with the nonlinear numerical solution for several initial amplitudes.
- Computes the nonlinear pendulum's period as a function of amplitude, and its relative error with respect to the linear model.
- Illustrates the **phase drift** between the two models.
- Simulates **energy dissipation** in a damped pendulum and compares it with the linear prediction.
- Generates the figures as PDF and PNG, and prints a ready-to-use LaTeX table.

## Mathematical model

Damped pendulum equation:

```
θ'' + c·θ' + ω²·sin(θ) = 0,      ω = √(g/L)
```

It is rewritten as a first-order system, `Y' = F(Y)` with `Y = (θ, v)`, and integrated with the RK4 scheme (local truncation error O(h⁵), global error O(h⁴)):

```
k₁ = h·F(Yₙ)
k₂ = h·F(Yₙ + k₁/2)
k₃ = h·F(Yₙ + k₂/2)
k₄ = h·F(Yₙ + k₃)
Yₙ₊₁ = Yₙ + (k₁ + 2k₂ + 2k₃ + k₄) / 6
```

**The period** is extracted from the numerical solution by detecting where the velocity `v` crosses zero (with linear interpolation between steps); the full period is twice this half-interval.

**Mechanical energy** (per unit mass and length²) is `E = v²/2 + ω²(1 − cos θ)`. For the damped pendulum, `dE/dt = −c·v² ≤ 0`, and in the linear regime the energy decays exponentially with time constant `τ = 1/c`.

## Simulation parameters

| Parameter | Value |
|---|---|
| Pendulum length `L` | 1 m |
| Gravitational acceleration `g` | 9.81 m/s² |
| Angular frequency `ω` | ≈ 3.132 rad/s |
| Linear period `T₀` | ≈ 2.006 s |
| Time step (simulations) | `h = 0.005 s` |
| Initial conditions | `θ(0) = θ₀`, `θ'(0) = 0` (released from rest) |
| Compared amplitudes | 10°, 60°, 90° |
| Damping coefficients | `c` = 0.2 / 0.5 / 1.0 s⁻¹ (θ₀ = 30°) |

## Results

> Note: the text inside the figures (titles, axis labels, legends) is in Romanian, as in the original thesis.

### 1. Linear model vs. RK4

![Linear vs RK4 comparison](fig1_comparatie.png)

At small amplitudes the two curves practically coincide. As the amplitude grows, the nonlinear pendulum's period becomes longer, and after a few periods the curves end up in antiphase (phase drift).

| Initial amplitude `θ₀` | Period error relative to the linear model |
|---|---|
| 10° | 0.19 % |
| 60° | 7.32 % |
| 90° | 18.03 % |

### 2. Period ratio `T(θ₀)/T₀` as a function of amplitude

![Period ratio T/T0](fig2_perioada.png)

### 3. Energy dissipation in the damped pendulum

![Energy dissipation](fig3_disipare.png)

Energy decays exponentially with time constant `τ = 1/c`: the larger `c` is, the faster the pendulum comes to rest.

### Conclusions

1. **Validity range of the linearization:** the approximation `sin θ ≈ θ` introduces period errors below 1 % only for `θ₀ < 20°`.
2. **Period lengthening:** at `θ₀ = 90°`, the real period is about 18 % longer than the linear prediction.
3. **Phase drift:** phase errors accumulate over time even when the amplitudes coincide.
4. **Energy dissipation:** the coefficient `c` directly determines the stopping time through `τ = 1/c`.

## Project structure

```
.
├── main.py                  # simulation, figure generation and LaTeX table
├── fig1_comparatie.png/pdf  # linear vs. RK4 comparison
├── fig2_perioada.png/pdf    # period ratio T(θ₀)/T₀ vs. amplitude
├── fig3_disipare.png/pdf    # energy dissipation (damped pendulum)
├── requirements.txt         # dependencies
└── README.md
```

## Installation

Requirements: **Python 3.9+**.

```bash
# 1. Clone the repository
git clone https://github.com/<username>/<repo>.git
cd <repo>

# 2. (Recommended) create a virtual environment
python -m venv .venv

# Windows
.venv\Scripts\activate
# Linux / macOS
source .venv/bin/activate

# 3. Install the dependencies
pip install -r requirements.txt
```

## Usage

```bash
python main.py
```

The script:

1. prints the simulation parameters (`L`, `g`, `ω`, `T₀`);
2. generates `fig1_comparatie`, `fig2_perioada` and `fig3_disipare` (PDF + PNG) in the current folder;
3. prints the table of RK4-computed periods to the console, together with the corresponding LaTeX code.

> Computing periods over a dense set of amplitudes uses small time steps, so a run may take a little while.

The parameters (`L`, `g`, the step `h`, the amplitudes, the damping coefficients) can be changed directly in `main.py`.

## Built with

- [Python](https://www.python.org/)
- [NumPy](https://numpy.org/) — numerical computation
- [Matplotlib](https://matplotlib.org/) — plotting

## Author

<!-- Fill in: First Last — GitHub / LinkedIn / email link -->

## License

<!-- Choose a license (e.g. MIT) and add a LICENSE file -->
