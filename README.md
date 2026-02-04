# Two-Locus Plant Mating Model (M/F/Fp/Fs/z)

This app simulates a diploid, hermaphroditic plant population with two tightly linked loci: an M locus (alleles M/m) and an F locus (alleles F/f/z/Fp/Fs). `Fp` is a rare mutation from `F` that is resistant to `z`-mediated reduction. `Fs` is an `F`-like allele with reduced efficacy (see selection rules) and is created only in diploid genomes that carry `z`. By default the loci are effectively nonrecombining, but you can optionally add a rare recombination rate between them.

## Mating behavior

Each generation is produced by random mating among all diploid genotypes built from the haplotypes. Plants are hermaphrodites, so any individual contributes both ovules and pollen. Fertilization is modeled by drawing a maternal genotype (for ovules) and then drawing a pollen haplotype from the population **with selection acting on pollen** depending on the maternal genotype.

### Pollen-side selection rules

Selection only changes **pollen success**, not ovule production.

- **m pollen penalty (`s_m`)**
  - If the maternal plant carries at least one **F**, **Fp**, or **Fs** allele, **m-bearing pollen** is disadvantaged by `(1 - s_m)`.
  - If the maternal genotype carries **z** and **F/Fp** (e.g., **Fz**, **Fpz**), the penalty is reduced to `s_m * (1 - z_reduction)` for that maternal genotype.
  - If the maternal genotype has **Fs** but **no F/Fp** (e.g., **FsFs** or **Fsz**), the penalty is reduced to `s_m * (1 - z_reduction)` (Fs is weaker unless dominated by F/Fp).
  - If the maternal plant has no F-like allele, there is **no** `s_m` penalty.

- **M pollen penalty (`s_M`)**
  - If the maternal plant is **ff**, **M-bearing pollen** is disadvantaged by `(1 - s_M)`.
  - If the maternal plant is not ff (i.e., has F or z), there is **no** `s_M` penalty.

These rules are applied **each time a maternal genotype is sampled**, so the pollen weights are conditional on the maternal genotype.

### Mutation

After mating and zygote formation, haplotypes mutate with one-way rates:

- `u_m`: M → m
- `u_f`: F → f
- `u_z`: f → z
- `u_fp`: F → Fp
- `u_fs`: F → Fs **only in diploid genomes that contain z** (each F copy in such a genome converts with probability `u_fs`)

No reverse mutation is modeled. Mutation is applied to haplotypes (not genotypes) and the haplotypes are renormalized afterward so they sum to 1.

### Recombination

After mutation, a recombination step moves haplotype frequencies toward linkage equilibrium at rate `r` between the M and F loci. Setting `r = 0` recovers the strictly linked case.

### Stopping rule

The simulation stops when both loci are effectively monomorphic (within tolerance) or when `max_generations` is reached.

## Mutation-only trajectories

The plot can optionally show dashed grey curves for how the overall **M** and **F** allele frequencies would decline **from mutation alone**, starting from the initial frequencies:

- `M(t) = M0 * (1 - u_m)^t`
- `F(t) = F0 * (1 - u_f - u_fp)^t`

`F0` is the initial frequency of **F-like** alleles (F and Fs, not `Fp` and not `z`). These are a visual reference and are not used in the simulation itself.

## Run from the command line

From the project directory:

```sh
R -q -e "shiny::runApp('.')"
```

This launches the Shiny app in your default browser.
