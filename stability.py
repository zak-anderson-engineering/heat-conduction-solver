import numpy as np
import matplotlib.pyplot as plt
from solver import solve_ftcs, analytical

if __name__ == "__main__":
    nx, t_end = 41, 100.0
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    for ax, r in zip(axes, [0.45, 0.55]):
        x, T = solve_ftcs(nx, t_end, r=r)
        print(f"r = {r}: largest temperature in solution = {np.max(np.abs(T)):.3e} C")

        ax.plot(x * 100, analytical(x, t_end), label="Analytical")
        ax.plot(x * 100, T, "o-", markersize=3, label=f"FTCS, r = {r}")
        verdict = "stable" if r <= 0.5 else "UNSTABLE"
        ax.set_title(f"r = {r}: {verdict}")
        ax.set_xlabel("Position along rod (cm)")
        ax.set_ylabel("Temperature above ends (°C)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"FTCS stability limit: r = $\\alpha \\Delta t / \\Delta x^2$ must be ≤ 0.5  (t = {t_end:.0f} s)")
    fig.tight_layout()
    fig.savefig("figures/stability.png", dpi=150)
    plt.show()