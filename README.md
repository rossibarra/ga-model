# Two-Locus Plant Mating Model (M/F/z)

This app simulates a diploid, hermaphroditic plant population with two tightly linked loci: an M locus (alleles M/m) and an F locus (alleles F/f/z). Tightly linked means **no recombination**, so haplotypes are inherited intact.

## Mating behavior

Each generation is produced by random mating among all diploid genotypes built from the haplotypes. Plants are hermaphrodites, so any individual contributes both ovules and pollen. Fertilization is modeled by drawing a maternal genotype (for ovules) and then drawing a pollen haplotype from the population **with selection acting on pollen** depending on the maternal genotype.

### Pollen-side selection rules

Selection only changes **pollen success**, not ovule production.

- **m pollen penalty (`s_m`)**
  - If the maternal plant carries at least one **F** allele (genotype includes F), **m-bearing pollen** is disadvantaged by `(1 - s_m)`.
  - If the maternal genotype is **Fz**, the penalty is reduced to `s_m * (1 - z_reduction)`.
  - If the maternal plant has no F allele, there is **no** `s_m` penalty.

- **M pollen penalty (`s_M`)**
  - If the maternal plant is **ff**, **M-bearing pollen** is disadvantaged by `(1 - s_M)`.
  - If the maternal plant is not ff (i.e., has F or z), there is **no** `s_M` penalty.

These rules are applied **each time a maternal genotype is sampled**, so the pollen weights are conditional on the maternal genotype.

### Mutation

After mating and zygote formation, haplotypes mutate with one-way rates:

- `u_m`: M → m
- `u_f`: F → f
- `u_z`: F → z

No reverse mutation is modeled. Mutation is applied to haplotypes (not genotypes) and the haplotypes are renormalized afterward so they sum to 1.

### Stopping rule

The simulation stops when both loci are effectively monomorphic (within tolerance) or when `max_generations` is reached.

## Mutation-only trajectories

The plot can optionally show dashed grey curves for how the overall **M** and **F** allele frequencies would decline **from mutation alone**, starting from the initial frequencies:

- `M(t) = M0 * (1 - u_m)^t`
- `F(t) = F0 * (1 - u_f - u_z)^t`

These are a visual reference and are not used in the simulation itself.

## Run from the command line

From the project directory:

```sh
R -q -e "shiny::runApp('.')"
```

This launches the Shiny app in your default browser.
