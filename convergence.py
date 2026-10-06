import numpy as np
import matplotlib.pyplot as plt
from solver import solve_ftcs, analytical, L

if __name__ == "__main__":
    t_end = 100.0
    nxs = [11, 21, 41, 81, 161, 321]     # each grid is twice as fine as the last

    dxs, errors = [], []
    for nx in nxs:
        x, T_num = solve_ftcs(nx, t_end, r=0.4)
        err = np.max(np.abs(T_num - analytical(x, t_end)))
        dxs.append(L / (nx - 1))
        errors.append(err)
    dxs = np.array(dxs)
    errors = np.array(errors)

    print("  nx    dx (mm)     max error (C)   ratio to previous")
    for i, (nx, dx, e) in enumerate(zip(nxs, dxs, errors)):
        ratio = errors[i - 1] / e if i > 0 else float("nan")
        print(f"{nx:4d}   {dx*1000:7.3f}     {e:.3e}        {ratio:.2f}")

    # Slope of log(error) vs log(dx) gives the order of accuracy
    slope, _ = np.polyfit(np.log(dxs), np.log(errors), 1)
    print(f"\nObserved order of accuracy: {slope:.2f} (expected 2)")

    ref = errors[0] * (dxs / dxs[0]) ** 2
    plt.loglog(dxs * 1000, errors, "o", markersize=6, label="Max error vs analytical")
    plt.loglog(dxs * 1000, ref, "--", label="Second-order reference ($\\propto \\Delta x^2$)")
    plt.xlabel("Grid spacing $\\Delta x$ (mm)")
    plt.ylabel("Max error (°C)")
    plt.title(f"Grid convergence at t = {t_end:.0f} s (observed order = {slope:.2f})")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.savefig("figures/grid_convergence.png", dpi=150)
    plt.show()