# Optimal BGP Route Reflector Placement in ISP Networks

## CPSC 4110/5110 — Facility Location Optimization — Course Project

---

## Overview

This project applies **facility location optimization** methods to the real-world networking problem of **BGP Route Reflector (RR) placement**. We formulate the problem as a **p-median problem** on real ISP network topologies and solve it using three methods:

1. **Integer Linear Program (ILP)** — exact solution via JuMP + HiGHS
2. **Greedy Heuristic** — fast forward-selection
3. **Local Search (Swap Heuristic)** — iterative improvement with random restarts

The project uses **real ISP topologies** from the Internet Topology Zoo dataset (260+ networks) and evaluates the methods across 5 comprehensive experiments.

---

## Quick Start

### 1. Run Full Experiments

```bash
julia run_experiments.jl
```

This runs all 5 experiments. Results are saved as CSV files in `results/`.

**Expected runtime:** 10–60 minutes depending on hardware (the ILP on larger topologies and local search on Kdl take the most time).

### 2. Generate Plots

```bash
julia run_visualizations.jl
```

Generates PNG plots in `results/plots/` for the report.

---

## Problem Description

### The BGP Route Reflector Problem

In large ISP networks, routers within the same Autonomous System (AS) use iBGP. A full mesh of iBGP sessions requires N(N-1)/2 connections — impractical at scale. **Route Reflectors (RRs)** reduce this: a few routers act as RRs, and all others (clients) peer only with the RR.

**The catch:** When an RR selects routes, it uses its own location to compute the best exit point (hot-potato routing via IGP metrics). This can force clients to use suboptimal exits — a client near LA might be forced to route through New York because the RR in New York selected that exit.

This is documented in **RFC 9107 (BGP Optimal Route Reflection)**.

### Mapping to Facility Location

| Concept | Networking Version |
|---|---|
| Facility | Route Reflector |
| Client / demand point | iBGP client router |
| Candidate location | Any router in the network |
| Cost | Routing suboptimality penalty |
| Decision | Binary: place RR at node i? |
| Problem type | **p-median** on a network graph |

---

## Mathematical Formulation

### Decision Variables
- `y_j ∈ {0,1}` — 1 if node j is selected as RR
- `x_{ij} ∈ {0,1}` — 1 if client i is assigned to RR j

### Objective
```
Minimize  Σ_i Σ_j  cost(i,j) × x_{ij}
```

Where `cost(i,j)` measures the routing penalty client i experiences when assigned to RR at j: averaged over all exit-point pairs, how much extra distance does client i travel because the RR at j picks an exit optimal for j rather than for i?

### Constraints
1. **Exactly k RRs:** `Σ_j y_j = k`
2. **Each client assigned once:** `Σ_j x_{ij} = 1  ∀i`
3. **Assign only to open RRs:** `x_{ij} ≤ y_j  ∀i,j`

---

## Solution Methods

### Method 1: ILP (Exact)
- Formulated using **JuMP** with the **HiGHS** solver (free, no license needed)
- Gives provably optimal solution
- Becomes slow on large networks (NP-hard), demonstrating the need for heuristics

### Method 2: Greedy Heuristic
- Forward selection: at each step, add the node reducing total cost the most
- Time: O(k × n²)
- Fast, typically within 5–15% of optimal

### Method 3: Local Search (Swap)
- Start with k random RRs
- Try all single-swap improvements; accept if objective decreases
- Multiple random restarts (best of 5)
- Good balance of speed and quality

### Baselines
- **Betweenness Centrality** — place RRs at graph-central nodes
- **Highest Degree** — place at most-connected nodes
- **Random** — best of 20 random placements

---

## Experiments

### Experiment 1: Optimality Gap
Compare ILP vs Greedy vs Local Search on small/medium topologies (k = 1,2,3,5).

### Experiment 2: Scaling Behavior
Run all methods on progressively larger topologies (11 → 754 nodes). Plot runtime and quality vs size.

### Experiment 3: Impact of k
Fix Cogentco (197 nodes), vary k from 1 to 20. Show diminishing returns.

### Experiment 4: Network Structures
Compare across US, European, and global topologies. Analyze how density, degree distribution, and structure affect optimal placement.

### Experiment 5: Optimized vs Naive
Quantify improvement of greedy/local search over betweenness, degree, and random baselines.

---

## Data

All data comes from the **Internet Topology Zoo** (Knight et al., IEEE JSAC 2011):
- URL: https://topology-zoo.org/dataset.html
- 260+ real ISP topologies in GML format
- Nodes have geographic coordinates (lat/lon)
- Edge weights computed via Haversine distance (proxy for OSPF/IGP costs)
- Exit points: degree-1 nodes (stub routers) supplemented to ~20% of nodes

---

## Technology Stack

- **Language:** Julia 1.9+
- **Optimization:** JuMP + HiGHS (free MIP solver)
- **Graphs:** Graphs.jl
- **Visualization:** Plots.jl + Colors.jl
- **Data output:** CSV.jl + DataFrames.jl

---

## References

1. E. Rosen, Y. Rekhter. BGP Route Reflection: An Alternative to Full Mesh Internal BGP (RFC 4456). IETF 2006.
2. R. Raszuk, B. Decraene, C. Cassar, K. Wang. BGP Optimal Route Reflection (RFC 9107). IETF, 2021.
3. D. Walton, A. Retana, E. Chen, J. Scudder. Advertisement of Multiple Paths in BGP (RFC 7911). IETF, 2016. 
4. S. Knight, H. Nguyen, N. Falkner, R. Bowden, M. Roughan. The Internet Topology Zoo. IEEE JSAC, 29(9), 2011. 
5. Internet Topology Zoo dataset: https://topology-zoo.org/dataset.html
6. Y. Breitbart, M. Garofalakis, A. Gupta, A. Kumar, R. Rastogi. On Configuring BGP Route Reflectors. 
IEEE COMSWARE, 2007. 
7. R. Benkoczi. Introduction to Facility Location Problems. CPSC 4110/5110 Lecture Notes, University of Lethbridge, 2026. 
8. R. Benkoczi. Node Search Algorithms for the p-Center. CPSC 4110/5110 Course Notes, 2026. 
9. G. Cornuejols, M. Fisher, G. Nemhauser. Facility Location. Integer Programming Tutorial, 1995. 
10. JuMP.jl: Modeling language for mathematical optimization in Julia. https://jump.dev/ 
11. HiGHS: High-performance open-source linear programming solver. https://highs.dev/ 
---
