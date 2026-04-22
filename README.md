# Fuzzy-Sphere Study of the Cubic CFT

This repository contains the Exact Diagonalisation (ED) and Density Matrix Renormalisation Group (DMRG) code used in our study of the **cubic fixed point**, regularised on the fuzzy sphere, as described in

> **Quantum Rotors on the Fuzzy Sphere and the Cubic CFT**  
> Andreas Stergiou (2026)  
> [arXiv:2604.xxxxx](https://arxiv.org/abs/2604.xxxxx)

## Physics Background

### The Fuzzy-Sphere Regularisation

The fuzzy sphere is a non-perturbative regularisation of 2+1D quantum field theories that preserves the full rotational SO(3) symmetry of the continuum. The starting point is a quantum rotor model: each of the `nm` sites on the sphere hosts a linear rigid rotor, whose orientation is specified by a unit vector $\hat{\mathbf{n}}_i$. Each rotor is truncated to its lowest two angular momentum shells, $\ell = 0$ (singlet) and $\ell = 1$ (triplet), giving four states per site. These four states are mapped onto four fermionic flavours, and the rotors are placed on a fuzzy sphere by projecting to the lowest Landau level of fermions in a magnetic monopole field. This gives `no = nm × 4` total fermionic modes.

When tuned to a critical point, the low-energy spectrum of the finite sphere directly encodes the **operator spectrum** of the corresponding (2+1)-dimensional CFT through the state-operator correspondence.

### The Hamiltonian

The full Hamiltonian is

$$\hat{H} = \hat{H}_{\mathrm{Hub}} + \hat{H}_{\mathrm{Heis}} + h\,\hat{H}_{\mathrm{trans}} + w\,\hat{H}_{\mathrm{cubic}}$$

| Term | Role |
|------|------|
| $\hat{H}_{\mathrm{Hub}}$ | Flavour-agnostic density-density repulsion (Haldane pseudopotentials `[6.5, 1.0]`); enforces single occupancy per orbital |
| $\hat{H}_{\mathrm{Heis}}$ | Heisenberg-like interaction via the orientation operator $\hat{\mathbf{n}}$, which connects the singlet and triplet flavours; realises the rotor-rotor coupling $\hat{\mathbf{n}}_i \cdot \hat{\mathbf{n}}_j$ on the fuzzy sphere (coupling fixed at 1.4) |
| $h\,\hat{H}_{\mathrm{trans}}$ | Single-particle energy splitting between singlet and triplet flavours, corresponding to the kinetic term of the rotor model; `h` is tuned to reach the critical point |
| $\hat{H}_{\mathrm{cubic}}$ | Cubic anisotropy built from projectors $P^x, P^y, P^z$ onto the Cartesian triplet states $\|x\rangle, \|y\rangle, \|z\rangle$; breaks O(3) down to the cubic group; coupling `w` is the control parameter |

### Symmetry Sectors

The Hilbert space is decomposed by four symmetry quantum numbers:

| Label | Meaning |
|-------|---------|
| Z | Total parity (0 or 1) — product of parities of all flavour occupations |
| P2 | N₂ parity (0 or 1) — parity of the number of flavor-2 particles |
| R | Spatial reflection on the sphere (±1); ED only |
| X | Flavour exchange 1 ↔ 3 (±1); ED only |

ED uses all 16 combinations of these four. DMRG uses only the 4 (Z, P2) sectors (ITensors supports at most 4 quantum numbers).

### Operators and Observables

For each eigenstate we measure the energy E, angular momentum $L^2 = \ell(\ell+1)$, and SO(3) Casimir $C^2$. Six key operators are tracked:

| Operator | $\ell$ | $C^2$ | Degeneracy |
|----------|--------|-------|------------|
| $\phi$ | 0 | ≈ 2 | 3 |
| $X$ | 0 | ≈ 6 | 2 |
| $Z$ | 0 | ≈ 6 | 3 |
| $A_\mu$ | 1 | ≈ 2 | 3 |
| $T_{\mu\nu}$ | 2 | ≈ 0 | 1 |
| $\epsilon'$ | 0 | ≈ 20 | 1 |

Energies are normalised so that $E_\mathrm{vac} = 0$ and $E_\phi = 0.5189$.

### Finite-Size Scaling of h

The polarisation parameter `h` is optimised for each (nm, w). It follows the scaling:

$$h_{\mathrm{opt}}(n_m) = h_c + \frac{A}{n_m^{b/2}}$$

| w | h_c | A | b |
|---|-----|---|---|
| 0.2 | 14.985 | 12.69 | 2.77 |
| 0.4 | 14.957 | 15.40 | 2.94 |
| 0.6 | 14.955 | 24.69 | 3.39 |
| 0.8 | 14.855 | 15.40 | 2.87 |

## Repository Structure

```
├── ed/
│   └── cubic_spectrum_ed.jl     # Full ED: all 16 sectors, combined output
│
├── dmrg/
│   ├── cubic_spectrum_dmrg_build_mpo.jl  # Build and cache H, L², C² MPOs
│   ├── cubic_spectrum_dmrg.jl            # DMRG for one (Z, P2) sector
│   ├── cubic_spectrum_dmrg_combine.jl    # Combine 4 sector results
│   └── run_dmrg.sh                       # Convenience script: runs all three steps
```

## Installation

### Julia packages

```julia
using Pkg
Pkg.add(["FuzzifiED", "ITensors", "ITensorMPS"])
```

- [**FuzzifiED**](https://github.com/FuzzifiED/FuzzifiED.jl) — fuzzy-sphere Hamiltonian builder and ED routines
- [**ITensors.jl**](https://github.com/ITensor/ITensors.jl) / **ITensorMPS.jl** — quantum-number-conserving DMRG

## Running ED

The ED script loops over all 16 symmetry sectors, diagonalises each, and writes the combined spectrum to a single output file.

```bash
cd ed/
julia cubic_spectrum_ed.jl [nm] [w] [h]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `nm` | 8 | System size |
| `w`  | 0.6 | Cubic anisotropy coupling |
| `h`  | from fit | Polarisation parameter (auto-computed if omitted) |

**Example:**
```bash
julia cubic_spectrum_ed.jl 12 0.6
```

**Output:** `spectrum_nm12_w0p6.txt` — one row per degenerate multiplet, with columns `ScaledE | L2 | C2 | Degen | Irrep`.

> **Runtime:** nm=8 takes a few seconds; nm=10 a few minutes; nm=12 up to ~30 minutes on a modern laptop.

## Running DMRG

The simplest way is via the shell script, which runs all three steps automatically:

```bash
cd dmrg/
./run_dmrg.sh <nm> <w> <h> [n_states] [threads]
```

**Example:**
```bash
./run_dmrg.sh 12 0.6 15.323
```

This runs in sequence:
1. **Build MPO cache** — builds the Hamiltonian and observable (L², C²) MPOs and saves them to `mpo_cache_nm12_w0p6_h15p323.jls`. Skipped if the cache already exists.
2. **DMRG sectors** — runs `cubic_spectrum_dmrg.jl` for each of the 4 (Z, P2) sectors.
3. **Combine** — merges the sector results into `dmrg_combined_nm12_w0p6_h15p323.txt`.

You can also run the steps individually:

```bash
# Step 1
julia cubic_spectrum_dmrg_build_mpo.jl 12

# Step 2 (one per sector)
julia cubic_spectrum_dmrg.jl 0 0 8 12   # Z=0, P2=0
julia cubic_spectrum_dmrg.jl 0 1 8 12   # Z=0, P2=1
julia cubic_spectrum_dmrg.jl 1 0 8 12   # Z=1, P2=0
julia cubic_spectrum_dmrg.jl 1 1 8 12   # Z=1, P2=1

# Step 3
julia cubic_spectrum_dmrg_combine.jl 12 15.323 0.6
```

### Adjusting parameters

In `cubic_spectrum_dmrg_build_mpo.jl` and `cubic_spectrum_dmrg.jl`, the parameters `h` and `w` are set near the top of the file (lines ~28–29). Change these to match your target (nm, w) run.

### Key DMRG features

- **Checkpointing** — a checkpoint is saved after every sweep. If a run is interrupted, re-run with `--resume` to continue from where it stopped.
- **Convergence** — early stopping triggers when both the energy *and* the mid-chain von Neumann entropy change less than the tolerance for two consecutive sweeps.
- **Excited states** — computed sequentially using a penalty-weight projector (weight = 100).
- **Subspace re-diagonalisation** — after all states converge, H is re-diagonalised in the DMRG subspace to remove bias from the penalty projector.

### Bond dimensions by system size

| nm | Max bond dim | Sweeps |
|----|-------------|--------|
| ≤ 12 | 1500 | 20 |
| 14 | 2000 | 26 |
| 16 | 2500 | 30 |
| ≤ 20 | 3000 | 36 |
| 22 | 3500 | 40 |

> **Runtime note:** nm=12 with 8 states takes a few hours on a modern workstation. Larger nm values require proportionally more memory and time.


## Output File Format

**ED** — `spectrum_nm<nm>_w<w>.txt`:
```
# ScaledE  | L2      | C2      | Degen | Irrep
  0.00000  |  0.000  |  0.000  |     1 | A1g
  0.51890  |  0.000  |  2.000  |     3 | T1u
  ...
```

**DMRG combined** — `dmrg_combined_nm<nm>_w<w>_h<h>.txt`:
```
ScaledE | Energy | L^2 | C^2 | Z | P2 | State#
```


## Citation

If you use this code, please cite the paper 

```bibtex
@article{Stergiou:2026jbw,
    author = "Stergiou, Andreas",
    title = "{Quantum Rotors on the Fuzzy Sphere and the Cubic CFT}",
    eprint = "2604.xxxxx",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    month = "4",
    year = "2026"
}
```

and the libraries

- FuzzifiED: [https://github.com/FuzzifiED/FuzzifiED.jl](https://github.com/FuzzifiED/FuzzifiED.jl)
- ITensors.jl: M. Fishman, S. R. White, E. M. Stoudenmire, *SciPost Phys. Codebases* **4** (2022)


## License

Released under the [MIT License](https://opensource.org/licenses/MIT).
