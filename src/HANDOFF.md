# HybridDecomposition handoff

## Current layout

- `src/HybridDecomposition.m2`: core numerical decomposition and exact-recovery code, followed by the test function.
- `src/ExampleSystems.m2`: exact `QQ` translations of the nine systems in `assets/examples/*.jl`, plus one toy system.
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

## Examples and focused test

`ExampleSystems.m2` exports all systems through `exampleSystems` and `exampleSystemNames`.

The added toy system is `toy-three-points-and-circle`:

```text
((x^2+y^2-1) x(x-1)(x-2), (x^2+y^2-1)(y-3)).
```

Its variety is the circle `x^2+y^2=1` together with `(0,3)`, `(1,3)`, and `(2,3)`.

`testHybridDecompositionOnExamples()` currently tests only this toy system. It checks:

- dimensions `{0,1}`;
- total degree `3` in dimension zero;
- total degree `2` in dimension one;
- successful exact reconstruction of all three isolated points.

The numerical curve component is present in `(result#"ByDimension")#1`. Its exact `f` reconstruction normally fails with `"no stable full-degree integer relation for f"` because NAG's stored generic slice is complex rather than rational. The top-level hash-table printer abbreviates `ByDimension`, and the current `PrintResults` block is reached only after `f` is recovered, so the failed curve is not printed as a separate `dimension 1` line.

## Rational slicing status

The requested post-decomposition rational reslicing code was removed. The current implementation still calls ordinary NAG `regeneration`, so positive-dimensional exact recovery is not guaranteed.

NAG has no public `RationalSlices => true` option. A possible future change is a custom regeneration wrapper that initializes `numericalAffineSpace` with an invertible integer matrix and integer constant vector, then applies `hypersurfaceSection` successively. The final stored slices inherit those rational rows. This approach was discussed and prototyped interactively, but it is **not implemented in the repository**. NAG may still use temporary complex slices internally for homotopy and monodromy.

## Running and current footer caveat

The files now live in `src`, and `HybridDecomposition.m2` loads `ExampleSystems.m2` by a relative filename. Run Macaulay2 from `src`, for example:

```sh
cd src
M2 -q
```

Then at the Macaulay2 prompt:

```m2
load "HybridDecomposition.m2";
exampleTestResults = testHybridDecompositionOnExamples();
```

The file currently contains `end;` before a manual `clearAll`/self-load/test footer. When loaded interactively, `end;` prevents that footer from running, which is why the test is called manually. Consequently, the `M2 --script HybridDecomposition.m2` command shown in the README does not execute the footer in the current state. A future cleanup should choose one workflow and make the README and footer consistent.

Several translated Julia systems have high path counts and can take a long time to decompose.
