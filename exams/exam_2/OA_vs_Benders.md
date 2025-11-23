## 🔹 The Setup: MINLP Formulation

A general MINLP has the form:

$$
\begin{aligned}
\min_{x, y} \quad & f(x, y) \\
\text{s.t.} \quad & g(x, y) \le 0, \\
& x \in \mathbb{R}^n, \; y \in \{0,1\}^m.
\end{aligned}
$$

where:

* $x$: continuous (real) variables  
* $y$: discrete (integer or binary) variables  
* $f(\cdot), g(\cdot)$: nonlinear, possibly nonconvex functions  

Both OA and BD exploit **decomposition** to separate the continuous and discrete parts, but they do so differently.

---

## ⚙️ 1. Outer Approximation (OA)

### **Core Idea**

OA alternates between solving:

1. A **NLP subproblem** (continuous relaxation with fixed integers)  
2. A **MILP master problem** that linearizes the nonlinear functions  

The nonlinearities are approximated **from the outside** (hence "outer approximation") using linearizations (cuts) of the nonlinear constraints and objective around previously computed NLP solutions.

---

### **Algorithm Sketch**

1. **Fix integer variables** $y = y^k$ and solve the **NLP subproblem:**

   $$
   \min_{x} f(x, y^k) \quad \text{s.t.} \quad g(x, y^k) \le 0.
   $$

   This gives continuous solution $x^k$ and objective value $z^k$.

2. **Build a MILP master problem** using:

   * Linear approximations (first-order Taylor expansions) of $f(x, y)$ and $g(x, y)$ around all previous NLP points $(x^i, y^i)$  
   * Integer constraints on $y$

3. **Solve MILP master problem** to get a new integer vector $y^{k+1}$.

4. **Repeat** until the upper (NLP) and lower (MILP) bounds converge.

---

### **Key Properties**

* Works best for **convex MINLPs**, since linearizations are valid outer approximations.  
* Each NLP is solved to *local* or *global* optimality given fixed $y$.  
* The master MILP accumulates linear cuts → successively tighter approximation.

---

## ⚙️ 2. Benders Decomposition (BD)

### **Core Idea**

Benders splits the problem into:

* A **master problem** involving only **integer variables** $y$  
* A **subproblem** involving continuous variables $x$  

Instead of linearizing nonlinear functions, it adds **Benders cuts** based on the dual information of the subproblem.

---

### **Algorithm Sketch**

1. **Fix integer variables** $y = y^k$ and solve the **NLP subproblem:**

   $$
   \min_{x} f(x, y^k) \quad \text{s.t.} \quad g(x, y^k) \le 0.
   $$

   Compute optimal value $z^k$.

2. **Formulate Benders cuts**:

   * **Feasibility cuts** if subproblem is infeasible for given $y^k$  
   * **Optimality cuts** if subproblem is feasible, derived from dual multipliers:

     $$
     \theta \ge \alpha^k + \beta^{kT}(y - y^k)
     $$

     where $\theta$ approximates the subproblem’s cost.

3. **Solve master problem:**

   $$
   \min_{y, \theta} \quad \theta \\
   \text{s.t. Benders cuts, } y \in \{0,1\}^m.
   $$

4. **Repeat** until convergence.

---

### **Key Properties**

* Focuses on *dual-based* cuts (not primal linearizations).  
* Suited for **separable** or **convex problems** where dual multipliers are available.  
* More scalable when the continuous subproblem is large but the integer part is small.  
* Extended to **Generalized Benders Decomposition (GBD)** for nonlinear cases.

---

## 🔸 Comparison Summary

| Feature                 | **Outer Approximation (OA)**                          | **Benders Decomposition (BD / GBD)**                    |
| ----------------------- | ----------------------------------------------------- | ------------------------------------------------------- |
| Master problem          | MILP with linearized constraints                      | MILP/MINLP with Benders cuts                            |
| Subproblem              | NLP with fixed $y$                                    | NLP with fixed $y$ (used for cut generation)            |
| Type of cuts            | Tangent plane linearizations (primal cuts)            | Dual-based cuts (feasibility/optimality)                |
| Requires convexity      | Yes (for guarantees)                                  | Yes (for convergence guarantees)                        |
| Typical use case        | Continuous nonlinearities are mild, easy to linearize | Large continuous subproblem, dual structure exploitable |
| Computational focus     | Approximating nonlinear feasible region               | Approximating cost/value function of subproblem         |
| Dual multipliers needed | No                                                    | Yes                                                     |
| Implementation examples | DICOPT (GAMS), Bonmin (OA mode)                       | GBD (generalized Benders), AIMMS Benders solvers        |

---

## 🔹 Intuition

* **OA:** builds a *piecewise-linear outer shell* of the nonlinear feasible set in the $(x, y)$ space.  
* **Benders:** builds a *piecewise-linear lower bound* on the *optimal value function* $\phi(y) = \min_x f(x,y)$.

---

## 🔹 In Practice

* **Outer Approximation** tends to outperform Benders on “balanced” MINLPs (moderate continuous and integer parts).  
* **Benders Decomposition** is better when the integer dimension is small but continuous subproblems are large and separable.
