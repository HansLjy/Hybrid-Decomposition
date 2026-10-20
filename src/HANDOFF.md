# HybridDecomposition handoff

## Current layout

- `src/HybridDecomposition.m2`: core numerical decomposition and exact-recovery code, followed by the test function.
- `src/ExampleSystems.m2`: exact `QQ` translations of the nine systems in `assets/examples/*.jl`, plus one toy system.
- `src/ToyExamples.m2`: small exactly solved rational slices used by the direct test loop.
- `src/README.md`: algorithm, output structure, LLL construction, multiplicity handling, and test notes.

The code is deliberately a plain Macaulay2 script, not a package.

## Implemented computation

`hybridDecomposition(I)` expects an ideal in a polynomial ring over `QQ`.

1. It calls NAG `regeneration` and, by default, asks for singular output.
2. Witness sets are separated by dimension and then set-theoretically by `decompose`.
3. For each component it stores the original NAG hyperplanes, witness set, witness points, and multiplicity diagnostics.
4. It chooses an integer linear form `t` separating the witness points.
5. It reconstructs an integer polynomial `f` whose roots are the projected values `t(alpha)`, with derivative equations included according to clustered numerical multiplicities.
6. It computes `fSquareFree = f // gcd(f,f')` and reconstructs each coordinate as

   `alpha_i = (h_i(t(alpha))/d_i) / fSquareFree'(t(alpha))`,

   with a separate integer denominator `d_i` and `degree(h_i) < degree(fSquareFree)`.

Results are returned in a hash table with `"Dimensions"` and `"ByDimension"`. The value `H#"ByDimension"#d` is a list because several irreducible components can have the same dimension.

## Where LLL is used

`relationsAtScale` constructs the column lattice

```text
[ I ; round(S Re(M)) ; round(S Im(M)) ]
```

and calls `LLL lattice`. The top block of each reduced vector is treated as a candidate integer relation. Candidates are filtered by height and relative residual.

`integerRelations` runs this at `S` and `S/100` and accepts only primitive relations found at both scales. This is a numerical stability heuristic, not a proof.

The same integer-kernel solver is used:

- once on the derivative/multiplicity matrix for all coefficients of `f` together;
- once per coordinate on the joint vector `(h_0,...,h_(s-1),d_i)`.

This is not coefficient-wise reconstruction.

## Multiplicity handling

The sliced system is solved numerically and clustered with `solutionsWithMultiplicity`. Cluster multiplicity is used at singular witness points. Regular witness points are assigned multiplicity one because duplicate total-degree paths are not necessarily local multiplicity.

## Examples and direct tests

`ExampleSystems.m2` exports the larger translated systems through
`exampleSystems` and `exampleSystemNames`.

`ToyExamples.m2` separately exports five directly constructed witness sets.
Every entry records `n`, `d`, an ideal over `QQ`, exactly `d` rational
hyperplanes, and all known slice solutions. The examples cover several
ambient dimensions below five and include irrational solutions involving
`sqrt(-2)` and `sqrt(3)`.

`testToyExamples()` loops over the full catalog. It verifies dimensions,
hyperplane counts, zero-dimensional slice degrees, residuals of all known
solutions, and successful exact reconstruction by `componentRepresentation`.

`testSys2449EquidimensionalComponents()` tests the dimension-eight,
degree-four top-dimensional part of `Sys2449`. It intersects the system with
eight seeded random rational hyperplanes, numerically groups the four slice
points into irreducible components of degrees one and three, reconstructs an
RUR for each group, and checks the RUR-derived solutions against both the
witness points and the original equations. `testExamples()` runs this test
after the toy catalog.

## Rational slicing status

The requested post-decomposition rational reslicing code was removed. The current implementation still calls ordinary NAG `regeneration`, so positive-dimensional exact recovery is not guaranteed.

NAG has no public `RationalSlices => true` option. A possible future change is a custom regeneration wrapper that initializes `numericalAffineSpace` with an invertible integer matrix and integer constant vector, then applies `hypersurfaceSection` successively. The final stored slices inherit those rational rows. This approach was discussed and prototyped interactively, but it is **not implemented in the repository**. NAG may still use temporary complex slices internally for homotopy and monodromy.

## Running and current footer caveat

The files now live in `src`, and `HybridDecomposition.m2` loads its support
files by relative filename. Run Macaulay2 from `src`, for example:

```sh
cd src
M2 -q
```

Then at the Macaulay2 prompt:

```m2
load "HybridDecomposition.m2";
testResults = testExamples();
```

Several translated Julia systems have high path counts and can take a long time to decompose.
