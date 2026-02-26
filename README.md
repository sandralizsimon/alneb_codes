# ALNEB: Active Learning Nudged Elastic Band

A MATLAB implementation of the Active Learning Nudged Elastic Band (ALNEB) method for finding minimum energy paths (MEPs) and transition states in atomistic systems. The method combines Gaussian Process Regression (GPR) with the Nudged Elastic Band (NEB) technique to efficiently locate transition states with minimal true energy/force evaluations.

## Method Overview

ALNEB operates in three phases:

### Phase 1: Exploration
Builds an initial GPR surrogate of the potential energy surface by iteratively:
1. Training a GPR model on known energy and force data
2. Running NEB optimization on the surrogate surface to find an approximate MEP
3. Selecting new training points based on maximum predicted uncertainty (covariance)
4. Evaluating true energies and forces at selected configurations
5. Repeating until the MEP movement between iterations falls below a convergence threshold

### Phase 2: Exploitation
Refines the MEP near the transition state region by:
1. Sampling at the highest predicted energy point along the MEP
2. Re-training the GPR and re-optimizing the NEB path
3. Converging until the highest energy image is accurately resolved
4. Performing a bracketing check to confirm the transition state candidate is a true energy maximum

### Phase 3: Renunciation
Performs precise transition state optimization using Climbing Image NEB (CI-NEB):
1. Extracts a local hypersphere of training data around the transition state guess
2. Trains a local GPR model on this subset
3. Runs CI-NEB on the surrogate surface to refine the climbing image
4. Evaluates true forces at the predicted TS geometry
5. Iterates until the maximum force component at the TS falls below 0.02 eV/A

## Requirements

- **MATLAB** (R2020a or later recommended)
- **Python 3** with the following packages:
  - [ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/)
- **Optimization Toolbox** (for `fmincon`)
- **Parallel Computing Toolbox** (optional, for faster hyperparameter optimization)

## File Structure

### Main Scripts
| File | Description |
|------|-------------|
| `exploration.m` | Phase 1 - Exploration loop: builds the initial GPR surrogate and finds an approximate MEP |
| `exploitation.m` | Phase 2 - Exploitation loop: refines the MEP and identifies the transition state region |
| `renunciation.m` | Phase 3 - Renunciation loop: precise TS optimization using CI-NEB on a local GPR model |

### GPR Module
| File | Description |
|------|-------------|
| `GPR_alpha_pinv_bayes.m` | Bayesian hyperparameter optimization and GPR weight computation |
| `GPR_predict_pinv.m` | GPR prediction of energies, forces, and covariances at new configurations |
| `rbf_kernel_extended_star.m` | Extended RBF kernel supporting energy, force, and cross-covariance terms |
| `negLogLikelihood_pinv_bayesian.m` | Negative log-likelihood objective for hyperparameter optimization |
| `get_conveged_surface.m` | Recomputes GPR weights using previously optimized hyperparameters (skip re-optimization) |

### NEB Module
| File | Description |
|------|-------------|
| `GPR_neb_opt_pinv.m` | NEB optimization on the GPR surrogate surface using an AARE optimizer |
| `GPR_ci_neb_opt_pinv.m` | Climbing Image NEB optimization on the GPR surrogate surface |
| `neb_force.m` | Standard NEB force projection (perpendicular true force + parallel spring force) |
| `cineb_force.m` | Climbing Image NEB force projection |
| `improved_tan.m` | Improved tangent estimation for NEB images |
| `spring_parallel.m` | Parallel component of spring forces between NEB images |

### Sampling and Analysis
| File | Description |
|------|-------------|
| `get_max_covariance.m` | Selects the next training point based on maximum GPR uncertainty |
| `find_high_energy.m` | Locates the highest energy configuration along the MEP on a finer grid |
| `finer_predict.m` | GPR predictions on a finer interpolation grid between NEB images |
| `finerdist.m` | Generates a finer grid of points between NEB images |

### I/O and Utilities
| File | Description |
|------|-------------|
| `readPOSCAR.m` | Reads VASP POSCAR format files |
| `writePOSCAR.m` | Writes VASP POSCAR format files |
| `find_forces_emt.m` | Evaluates true energies and forces using the EMT potential via ASE |
| `get_emt.py` | Python script that calls ASE's EMT calculator |
| `interpolate.m` | Linear interpolation between two endpoint configurations |
| `compute_perturbed_points.m` | Generates perturbed configurations near the endpoint minima |
| `create_new_poscar.m` | Reorders atoms by mobility (free atoms first, fixed atoms last) |
| `pos_to_points.m` | Flattens a position matrix (N x 3) to a row vector |
| `points_to_pos.m` | Reshapes a row vector back to a position matrix (N x 3) |
| `distance.m` | Computes maximum atomic displacement between two NEB paths |
| `extract_index.m` | Returns the active DOF indices for GPR training |
| `fix_atoms.m` | Returns the indices of atoms held fixed during NEB |

### Input Files
Heptamer Island Transition A : DOF 189
| File | Description |
|------|-------------|
| `POSCAR_initial.1` | Reactant (initial state) geometry in VASP POSCAR format |
| `POSCAR_final.1` | Product (final state) geometry in VASP POSCAR format |

## Usage

1. Place the reactant and product geometries as `POSCAR_initial.1` and `POSCAR_final.1` in the working directory.

2. Configure system-specific parameters in `exploration.m`:
   - EMT potential parameters (`De`, `a`, `re`, `rc`) or replace `find_forces_emt.m` with your own energy/force calculator
   - GPR hyperparameter bounds and initial values
   - Dimensionality of the system (`dimension`)
   - Convergence thresholds (`dmax_limit`, `rmax_limit`, `dmin_limit`)
   - Number of NEB images (set in `interpolate.m`, default: 7)

3. Configure atom indexing in:
   - `extract_index.m` - set the active degrees of freedom
   - `fix_atoms.m` - set the indices of constrained atoms

4. Run the calculation:
   ```matlab
   exploration
   ```
   The code will automatically transition through all three phases.

## Output Files

| File | Description |
|------|-------------|
| `finaldata.txt` | All sampled atomic configurations |
| `finalvalues.txt` | Corresponding energies and forces |
| `array_out.mat` | Iteration-level log (convergence metrics, MEP paths, predicted energies) |
| `ts_data_final.mat` | Transition state optimization history (coordinates, gradients, MEP) |
| `maxforce_gpr.txt` | Maximum NEB force per optimizer step |
| `normforce_gpr.txt` | NEB force norm per optimizer step |
| `conf*/` | Individual configuration directories with POSCAR files and calculator outputs |

## Key Algorithmic Details

- **GPR Kernel**: Extended RBF kernel that jointly models energies and forces with their cross-covariances, enabling learning from both energy and gradient observations.
- **NEB Optimizer**: Adaptive Accelerated Relaxation Engine (AARE) combines Conjugate Gradient direction scheme (Fletcher-Reeves / Hestenes-Stiefel) with MD and adaptively accelerates based on the angle between successive gradients.
- **Active Learning**: New training points are selected either by maximum predicted covariance (uncertainty sampling) during exploration, or by maximum predicted energy (exploitation sampling) near convergence.
- **Convergence**: The MEP is considered converged when the maximum atomic displacement between successive iterations falls below `dmax_limit`. The transition state is converged when the maximum force component drops below 0.02 eV/A.
