# impact-beta

Numerical solver and parameter-space maps for the maximum spreading factor **β** of soft elastic (PAAm hydrogel) drops on rigid substrates, accompanying:

> Chowdhury, A., Mitra, S., & Mitra, S. K. *Bridging Liquid and Elastic Solid Impact Regimes Using Flexible Hydrogels.* Langmuir (2026). [doi:10.1021/acs.langmuir.6c02919](https://doi.org/10.1021/acs.langmuir.6c02919)

## What this repo does

The paper derives an energy-balance model for a spherical elastic drop (radius `r₀`, impact velocity `v₀`, shear modulus `G`, surface tension `γ`, density `ρ`) that deforms into a pancake of radius `βr₀` on impact. Balancing kinetic energy against neo-Hookean strain energy, surface energy, work of adhesion, and dissipation yields a **sixth-order polynomial in β** (Eq. 4):

```
β⁶ − [(1.5 + 0.5(1−λ)/El) / (1 − 0.75/(We·El))] · β⁴
    + [4 / (2·We·El − 1.5)] · β³
    + [0.5 / (1 − 0.75/(We·El))]  =  0
```

with

- `El = G/(ρv₀²)`  — elastic number
- `We = ρv₀²r₀/γ` — Weber number
- `λ = L/T`         — fraction of kinetic energy dissipated (empirical)

This notebook solves that polynomial numerically over a grid in (We, El), extracts positive real roots, and visualises the resulting β-landscape — including the two physically admissible roots (`β_min`, `β_max`) that appear in the intermediate-elasticity regime.

## Contents

| File | Description |
|---|---|
| `impact_beta.ipynb` | Main notebook: root-finder + all plots |
| `beta_max_map.csv`, `beta_min_map.csv` | β-solutions on the (We, El) grid |
| `We_vals.csv`, `El_vals.csv` | Grid axes |
| `min_beta_surface_We_El.png`, `max_beta_surface_We_El.png` | 3D surface plots |

## Method

For each `(We, El)` grid point:

1. Assemble the six polynomial coefficients from Eq. (4), skipping singular points where `We·El → 0.75` or `2·We·El → 1.5`.
2. Solve with `numpy.roots` (companion-matrix eigenvalues).
3. Retain roots with `|Im(β)| < 10⁻⁶·|Re(β)|` and `Re(β) > 0`.
4. Store `min` and `max` of the surviving roots.

Grid: `We ∈ [1, 500]`, `El ∈ [10⁻³, 10²]`, 500 × 500 log-spaced points. Default `λ = 0.5`.

## Reproducing the paper figures

The three impact velocities in the experiments (1, 2, 3 m/s) correspond to `We ≈ 30, 120, 270`. The notebook plots `β_max(El)` and `β_min(El)` slices at these values — these are the theoretical curves overlaid on the experimental data in Fig. 5(b,c) of the paper.

For the stiff-gel limit (`El ≫ 1`, `Ec ≫ 100`), the physical solution reduces to the cubic (Eq. 5), which the sixth-order roots recover as expected. For soft gels (`El < 1`), setting `λ ≈ 0.9–0.99` reproduces the purple/red dashed curves in Fig. 5.

## Usage

```bash
git clone <this repo>
cd impact-beta
pip install numpy matplotlib
jupyter notebook impact_beta.ipynb
```

Tune the top cell:

```python
LAMBDA = 0.5   # dissipation fraction, [0, 1)
N_WE   = 500   # We-axis resolution
N_EL   = 500   # El-axis resolution
```

Note: cell 5 contains a `google.colab` mount for saving CSVs to Drive — comment it out for local use.

## Physical interpretation

- **`El > 1`**: bulk elasticity dominates; `β_max = f(El)` alone (substrate-independent), and the neo-Hookean cubic matches experiment.
- **`El < 1`**: contact-line pinning and contact-foot ejection dissipate 90–99% of the kinetic energy; the polynomial reproduces the observed spread only with large `λ` and moderate `Ec`.
- The transition at `El ≈ 1` is where the impact force scaling also crosses over from the Wagner limit (`F* ≈ 3.24`) to the Hertzian power law (`F* ∼ El⁰·³⁸`).

## Citation

```bibtex
@article{Chowdhury2026Bridging,
  author  = {Chowdhury, Akash and Mitra, Surjyasish and Mitra, Sushanta K.},
  title   = {Bridging Liquid and Elastic Solid Impact Regimes Using Flexible Hydrogels},
  journal = {Langmuir},
  year    = {2026},
  doi     = {10.1021/acs.langmuir.6c02919}
}
```


