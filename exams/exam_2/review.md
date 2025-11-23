## **1. Augmented Lagrangian Method (ALM)**

**Intuition:**
Think of ALM as *Lagrange multipliers with a shock absorber.* The ordinary Lagrangian tries to enforce constraints through multipliers alone, but this can lead to slow convergence or oscillation. ALM adds a quadratic penalty that pulls iterates closer to feasibility **without using extremely large penalties**.

**Standard problem:**
$$
\min_x f(x) \quad \text{s.t. } h(x)=0
$$

**Ordinary Lagrangian:**
$$
\mathcal{L}(x,\lambda)= f(x)+\lambda^\top h(x)
$$

**Augmented Lagrangian:**
$$
\mathcal{L}_A(x,\lambda,\rho)=
f(x)+\lambda^\top h(x) + \frac{\rho}{2}|h(x)|^2
$$

**How it works intuitively:**

* The **quadratic term** penalizes constraint violation smoothly.
* The **dual variable update** enforces constraints like in Lagrange multipliers.
* The method alternates between minimizing in $x$ and updating multipliers:
  $$
  \lambda^{k+1} = \lambda^k + \rho h(x^{k+1})
  $$

**Why it's powerful:**

* Unlike penalty methods, you *don’t* need $\rho \to \infty$.
* Unlike pure dual ascent, the penalty term stabilizes convergence.
* Works well for nonconvex problems and forms the basis of **ADMM**.

---

## **2. Classic Integer Optimization Problems**

### **A. Knapsack Problem**

**Intuition:**
You have a bag with capacity $C$. Each item has value $v_i$ and weight $w_i$.
Pick items to maximize value without exceeding capacity.

**Formulation (0–1 knapsack):**
$$
\max_{x\in{0,1}^n} \sum_i v_i x_i \\
\text{s.t. } \sum_i w_i x_i \le C
$$

**Why it’s important:**

* Canonical NP-hard problem.
* Models resource allocation with limited budget.
* Basis of DP methods, branch & bound, and greedy approximations.

---

### **B. Traveling Salesman Problem (TSP)**

**Intuition:**
A salesman must visit every city once and return home using the shortest possible tour.

**Graph model:**
Nodes = cities, edges = distances.

**Formulation:**
Binary variable $x_{ij}=1$ if edge $i\to j$ is used.

$$
\min \sum_{i}\sum_{j} c_{ij} x_{ij}
$$

**Degree constraints:**
Every city has exactly one incoming and one outgoing edge:
$$
\sum_{j} x_{ij}=1,
\qquad
\sum_{j} x_{ji}=1
$$

**Subtour elimination constraints (SEC):**
$$
\sum_{i,j\in S} x_{ij} \le |S|-1,
\quad \forall S\subset{1,\ldots,n}
$$

**Importance:**

* NP-hard.
* Benchmark for branch-and-cut techniques.
* Highlights combinatorial explosion and need for relaxation/heuristics.

---

## **3. Outer Approximation (OA)**

**Intuition:**
OA is for **Mixed-Integer Nonlinear Programming (MINLP)** when the nonlinear part is convex. It alternates between:

1. solving a **NLP** for fixed integer variables,
2. building **linear cuts** (outer approximations) of the nonlinear functions,
3. solving a **MILP master problem** that uses those linearizations.

**Core idea:**
Approximate the nonlinear feasible region from the *outside* using linear cuts — the master MILP rules out infeasible or suboptimal integer choices.

**Algorithm steps (high level):**

1. Fix integer variables → solve NLP → get point $x^k$
2. Linearize constraints/objective at $x^k$:
   $$
   f(x) \ge f(x^k) + \nabla f(x^k)^\top (x - x^k)
   $$
3. Add these cuts to a growing MILP master problem.
4. Solve the master → pick new integer vector.
5. Repeat until bounds close.

**Why OA works:**

* NLP handles the hard nonlinear physics.
* MILP handles the discrete logic.
* Cuts gradually shrink the space of feasible discrete choices.

---

## **4. Benders Decomposition**

**Intuition:**
Benders separates a large optimization problem into:

* a **master problem** with complicating variables, and
* a sequence of smaller **subproblems** that enforce feasibility or compute an optimality cut.

**Typical form:**
$$
\min_x c^\top x + q(y) \\
\text{s.t. } Ax + By \ge b
$$

Where $y$ is “easy,” but the combination of $x$ $and $y$ is expensive.

**Key idea:**
Fix the master variables $x$, solve the subproblem (over $y$):

* If subproblem feasible → generate an **optimality cut**.
* If subproblem infeasible → generate a **feasibility cut** that excludes the current $x$.

The master accumulates cuts that tighten its feasible set.

**Mathematically:**
If dual of subproblem gives multipliers $\pi$, you add a cut:
$$
\theta \ge \pi^\top (b - A x)
$$

**Use cases:**

* Large-scale LP/MILP.
* Network design, logistics, energy systems.
* Problems with complicating constraints linking subsystems.

---

## **5. OA vs. Benders: Key Differences**

| Aspect         | Outer Approximation (OA)                           | Benders Decomposition                               |
| -------------- | -------------------------------------------------- | --------------------------------------------------- |
| Used for       | MINLP (convex in continuous vars)                  | LP/MILP with complicating constraints               |
| Master problem | MILP with linearizations of nonlinear constraints  | Simplified problem in master vars with Benders cuts |
| Subproblem     | NLP with fixed integer variables                   | LP/MILP over non-master vars                        |
| Cuts represent | Linearization of nonlinear functions               | Feasibility/optimality from dual subproblem         |
| Works when     | Nonlinearity is convex                             | Dual of subproblem easily available                 |
| Intuition      | “Approximate nonlinear region with tangent planes” | “Project out difficult variables using duality”     |

**High-level:**

* **OA**: linearizes convex nonlinearities.
* **Benders**: decomposes based on duality structure.
* Both alternate between *master problem* and *subproblem* but rely on different mathematics.