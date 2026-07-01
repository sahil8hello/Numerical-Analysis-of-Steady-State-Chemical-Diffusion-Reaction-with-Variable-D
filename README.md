# Steady-State Chemical Diffusion-Reaction with Variable D

This repository contains a MATLAB implementation that solves a steady-state mass balance ordinary differential equation (ODE) describing a chemical species $A$ diffusing into a liquid layer. The species simultaneously undergoes an irreversible, first-order chemical reaction that converts it to product $B$. 

The core highlight of this project is solving the system without making the simplifying assumption of a constant diffusion coefficient. Instead, it handles a **linearly varying diffusion coefficient $D(x)$** both analytically (via special mathematical functions) and numerically.

---

## 🔬 Problem Physics & Governing Equations

The system is described by the following 1D steady-state transport equation:

$$\frac{d}{dx} \left( D(x) \frac{dC_A}{dx} \right) - kC_A = 0, \quad 0 < x < L$$

### System Parameters
* **Domain Length ($L$):** $5.0 \text{ cm}$
* **Reaction Rate ($k$):** $0.1 \text{ /day}$ (First-order)
* **Variable Diffusion Coefficient ($D(x)$):** Varies linearly from $D_0 = 0.84 \text{ cm}^2\text{/day}$ at $x=0$ to $D_L = 1.09 \text{ cm}^2\text{/day}$ at $x=L$.
    $$D(x) = 0.84 + 0.05x$$

### Boundary Conditions
* **Constant Source ($x = 0$):** $C_A = 100 \text{ mol/L}$
* **No-Flux Impermeable Wall ($x = L$):** $\frac{dC_A}{dx} = 0$

---

## 🛠️ Methodologies Implemented

### 1. Exact Analytical Solution
By transforming variables into the spatial metric $z = \alpha\sqrt{D(x)}$, the variable-coefficient ODE is mapped directly to a standard **Modified Bessel Equation of Order Zero**:

$$\frac{d^2C_A}{dz^2} + \frac{1}{z} \frac{dC_A}{dz} - C_A = 0$$

The exact concentration profile is then mapped using Modified Bessel Functions of the First ($I_0$) and Second Kind ($K_0$):
$$C_A(x) = C_1 I_0(z(x)) + C_2 K_0(z(x))$$
*The integration constants $C_1$ and $C_2$ are explicitly solved using the system's boundary conditions via MATLAB's built-in `besseli` and `besselk` algorithms.*

### 2. Numerical Solution (Finite-Difference Method)
A second-order accurate central difference scheme [$O(\Delta x^2)$] is built by discretizing the derivative operator directly on half-steps to properly account for the varying spatial flux:

$$D_{i-\frac{1}{2}} C_{A,i-1} - \left(D_{i-\frac{1}{2}} + D_{i+\frac{1}{2}} + k\Delta x^2\right)C_{A,i} + D_{i+\frac{1}{2}} C_{A,i+1} = 0$$

* **Grid spacing ($\Delta x$):** $0.05 \text{ cm}$ (yielding $101$ discrete nodes)
* **Neumann Boundary ($x=L$):** Resolved with a second-order ghost-node formulation.

### 3. Natural Cubic Spline Interpolation
To capture continuity between the discrete $101$ spatial points, a continuous **Natural Cubic Spline** interpolation profile is constructed, explicitly enforcing that the second derivative of concentration vanishes at the outer physical boundaries ($\frac{d^2C_A}{dx^2} = 0$).

---

## 🚀 Getting Started

### Prerequisites
* MATLAB (R2021a or newer recommended)
* No external toolboxes required.

### Execution
1. Clone the repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/variable-diff-reaction-fdm.git](https://github.com/YOUR_USERNAME/variable-diff-reaction-fdm.git)