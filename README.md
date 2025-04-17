# Finite Differences Method for the 2D Wave Equation

This project implements the **Finite Differences Method (FDM)** to solve the 2D wave equation using both **explicit** and **implicit** numerical schemes in MATLAB/Octave. It also includes visualization of the simulation results, computation of the exact solution for comparison, and graphical analysis of the error.

---

## 📌 Features

- ✅ Explicit and implicit solvers for the 2D wave equation.
- 📊 Real-time 3D surface animations of wave propagation.
- 📉 Error visualization over time and mesh grid.
- 🔬 Support for custom initial and boundary conditions.
- 🧠 Includes exact analytical solution for benchmarking.
- 📈 Gauss-Seidel solver for the implicit method.

---

## 🧠 Method Explanation

This project solves the 2D wave equation of the form:

\[
\frac{\partial^2 u}{\partial t^2} = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) + f(x,y,t)
\]

using finite difference discretization over a rectangular domain \([0, l_j] \times [0, l_i]\) and applying:

- **Explicit Method:** A leapfrog-style method based on central differences in time and space.
- **Implicit Method:** A Gauss-Seidel iterative solver applied to the discretized system, more stable for large time steps.

The numerical solution is compared against the exact analytical solution, when provided, and error is computed and visualized accordingly.

---

## 📂 Project Structure

Each `.m` file serves a specific function:

| File                   | Purpose                                                        |
|------------------------|----------------------------------------------------------------|
| `fdm_explicit_animation.m` | Solves the 2D wave equation using the **explicit method**.     |
| `fdm_implicit_animation.m` | Solves the 2D wave equation using the **implicit method**.     |
| `exact.m`              | Computes the **exact analytical solution** if known.           |
| `error_graph.m`        | Generates a 3D surface plot of the **absolute error**.         |
| `gaussSeidel.m`        | Solves linear systems using the **Gauss-Seidel iteration**.    |

---

## 🧰 Requirements

- MATLAB **or** GNU Octave
- No additional toolboxes required
- Recommended: system with graphical support for real-time plotting

---

## 🚀 Installation

1. Clone this repository or download all `.m` files into a working directory:
    ```bash
    git clone https://github.com/arielpincayy/fdm-wave-2d.git
    cd fdm-wave-2d
    ```

2. Open MATLAB or Octave and set the current directory to this folder.

---

## 🔍 Usage

You can run either the explicit or implicit simulation by calling the corresponding function in MATLAB/Octave with the correct parameters.

### ✅ Example 1: Sinusoidal Initial Condition (Analytical)

```matlab
top = zeros(1,100);
bottom = zeros(1,100);
left = zeros(98,1);
right = zeros(98,1);

l_i = 3;
l_j = 3;

n = 100;
m = 100;

u_i = @(x,y,t) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
f = @(x,y,t) 0;

delta_t = 0.01;
t_end = 0.05;
c = 1;
c_w = 1;
z_range = [-1,1];

fdm_explicit_animation(top, bottom, left, right, c, u_i, f, l_i, l_j, delta_t, t_end, z_range);
fdm_implicit_animation(top, bottom, left, right, c_w, u_i, f, l_i, l_j, delta_t, t_end, z_range);

---

### ✅ Example 2: Polynomial Source Function

top = zeros(1,100);
bottom = zeros(1,100);
left = zeros(98,1);
right = zeros(98,1);

l_i = 3;
l_j = 3;

n = 100;
m = 100;

c = 1;
c_w = 1;

delta_t = 0.01;
t_end = 0.05;

f = @(x, y, t) (2 / c^2) * (x^2 - l_j * x) * (y^2 - l_i * y) ...
              - 2 * (y^2 - l_i * y) * t^2 ...
              - 2 * (x^2 - l_j * x) * t^2;

u_i = @(x, y, t) (x^2 - l_j * x) * (y^2 - l_i * y) * t^2;

z_range = [-1,3];

fdm_explicit_animation(top, bottom, left, right, c, u_i, f, l_i, l_j, delta_t, t_end, z_range);
fdm_implicit_animation(top, bottom, left, right, c_w, u_i, f, l_i, l_j, delta_t, t_end, z_range);

# 🧮 Gauss-Seidel Iteration (Implicit Solver)

The `gaussSeidel.m` function solves the linear system generated by the implicit discretization using iterative refinement:

```matlab
[x, iter] = gaussSeidel(A, b, tol, maxIter);
```
---
## Parameters

- `A`: Coefficient matrix (sparse or full)  
- `b`: Right-hand side vector  
- `tol`: Convergence tolerance (e.g., `1e-6`)  
- `maxIter`: Maximum number of iterations  

## Returns

- `x`: Solution vector  
- `iter`: Number of iterations performed  

---

## 📊 Output and Visualization

- 🌊 3D surface plots of `u(x, y, t)` at each timestep  
- 🧮 Error surface comparison between exact and approximate solution  
- 📈 Plot of error over time  
- 🔁 Plot of Gauss-Seidel iterations per timestep  

---

## ✍️ Author

Developed by arielpincayy. Contributions welcome!

