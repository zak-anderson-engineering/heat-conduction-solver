import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ALPHA = 1.2e-5   # thermal diffusivity, m^2/s (roughly steel)
L = 0.1          # rod length, m
T0 = 100.0       # peak starting temperature above the ends, deg C


def initial_condition(x):
    """Starting temperature: a sine-shaped hot spot, ends held at 0."""
    return T0 * np.sin(np.pi * x / L)


def analytical(x, t):
    """Exact solution for this case: the sine shape just decays exponentially."""
    return T0 * np.sin(np.pi * x / L) * np.exp(-ALPHA * np.pi**2 * t / L**2)


def solve_ftcs(nx, t_end, r=0.4):
    """Solve the 1D heat equation with FTCS (forward time, centred space).

    r = alpha * dt / dx^2 must be <= 0.5 or the solution becomes unstable.
    """
    x = np.linspace(0, L, nx)
    dx = x[1] - x[0]
    dt = r * dx**2 / ALPHA
    n_steps = max(1, int(np.ceil(t_end / dt)))
    dt = t_end / n_steps                  # adjust so we land exactly on t_end
    r = ALPHA * dt / dx**2

    T = initial_condition(x)
    for _ in range(n_steps):
        # each interior point is updated from its two neighbours; ends stay at 0
        T[1:-1] = T[1:-1] + r * (T[2:] - 2 * T[1:-1] + T[:-2])
    return x, T


if __name__ == "__main__":
    Path("figures").mkdir(exist_ok=True)
    nx = 41

    for t in [0, 50, 100, 200]:
        x, T_num = solve_ftcs(nx, t)
        T_exact = analytical(x, t)
        print(f"t = {t:3d} s   max error = {np.max(np.abs(T_num - T_exact)):.4f} C")
        line, = plt.plot(x * 100, T_exact, label=f"t = {t} s")
        plt.plot(x * 100, T_num, "o", markersize=3, color=line.get_color())

    plt.xlabel("Position along rod (cm)")
    plt.ylabel("Temperature above ends (°C)")
    plt.title(f"FTCS solver (points) vs analytical solution (lines), {nx} grid points")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("figures/solution_vs_analytical.png", dpi=150)
    plt.show()
