# HybridDecomposition

`HybridDecomposition.m2` is a plain Macaulay2 script that computes
component-separated numerical witness sets
and attempts to recover a rational univariate representation for each one.
The input must be an ideal in `QQ[x_1,...,x_n]`.

For each numerical irreducible component it returns:

- the numerical slicing hyperplanes and witness set returned by the
  decomposition;
- an integer separating linear form `t` and its projected witness values;
- an integer polynomial `f`, with clustered numerical multiplicities;
- `squarefree(f)`; and
- a separate rational polynomial `g_i=h_i/d_i` for every coordinate, where
  `degree(h_i) < degree(squarefree(f))` and
  `alpha_i = g_i(t(alpha))/squarefree(f)'(t(alpha))`.

## Use

```m2
load "HybridDecomposition.m2";
R = QQ[x,y];
I = ideal(y-x^2);
H = hybridDecomposition I;

-- There can be several components of the same dimension.
C = first ((H#"ByDimension")#1);
C#"Hyperplanes"
C#"WitnessPoints"
C#"f"
C#"SquareFreePart"
C#"NumeratorPolynomials"
C#"Denominators"
```

`H#"ByDimension"#d` is a list, so distinct numerical irreducible components
are never combined merely because they have the same dimension.

The program does not replace or move the slicing plane. Positive-dimensional
components are reconstructed directly from the numerical witness set and
hyperplanes produced by `regeneration` and `decompose`. Consequently, an
integer relation may legitimately be absent when a numerical slice is not
defined over `QQ`; in that case `Reconstructed` is `false` and the diagnostic
data are retained.

## How LLL is used

This is not coefficient-wise rational reconstruction. Every polynomial is
recovered by solving one approximate integer-kernel problem involving all of
its coefficients.

For `f=c_0+...+c_N T^N`, a projected value `beta_j` of numerical multiplicity
`m_j` supplies the simultaneous linear equations

```text
sum_r c_r (d^k/dT^k T^r)(beta_j) = 0,   0 <= k < m_j.
```

These rows form a complex matrix `M`. For a large integer scale `S`, the
columns of the integer lattice

```text
[ I ; round(S Re(M)) ; round(S Im(M)) ]
```

are LLL-reduced. A short lattice vector has top block `c` and bottom blocks
approximately `S Re(Mc)` and `S Im(Mc)`, so a short vector with small residual
provides the integer relation `Mc approximately 0`. The implementation
requires the same primitive relation to occur at both `S` and `S/100` before
accepting it.

After computing `f_red=squarefree(f)`, each coordinate is one more shared
integer-kernel problem. With `s=degree(f_red)`, its unknown vector is

```text
(h_0,...,h_(s-1),d_i)
```

and every distinct witness point supplies

```text
sum_r h_r beta_j^r - d_i alpha_(j,i) f_red'(beta_j) = 0.
```

Thus each coordinate gets its own numerator and denominator, but the
coefficients of that numerator and denominator are found together.

## Numerical multiplicity

The sliced polynomial system is solved numerically before reconstruction, and
nearby roots are clustered with `solutionsWithMultiplicity`. Cluster size is
used at witness points classified as singular. At a regular witness point the
multiplicity is fixed to one, because a total-degree homotopy can return
duplicate paths at a simple solution; those duplicates are not local
multiplicity. `ClusterTolerance` controls the clustering threshold and
defaults to `10^(-6)`.

The full diagnostic data are retained as `RawSliceSolutions`,
`ClusteredSliceSolutions`, `SingularWitnessPoints`, `FRelationData`, and the
per-coordinate `RelationData`.

Useful controls include:

```m2
hybridDecomposition(I,
    RelationScale => 10^60,
    RelationTolerance => 10^(-12),
    RelationMaxHeight => 10^30,
    RefinementBits => 300,
    RefinementTolerance => 10^(-70),
    ClusterTolerance => 10^(-6),
    SeparateComponents => true,
    PrintResults => true)
```

An exact result is refused if no full-degree relation is stable at both LLL
scales, if the squarefree degree does not equal the number of distinct witness
points, or if any coordinate relation fails.

## Direct witness-set tests and supplied examples

`ToyExamples.m2` contains five small ideals over $\mathbb{Q}$. Each entry
records its ambient variable count $n$, dimension $d$, ideal, exactly $d$
rational slicing hyperplanes, explicit known solutions, and the resulting
witness set. The examples cover
$(n,d)=(2,1),(3,1),(3,2),(4,2)$ and include slices whose solutions involve
$\sqrt{-2}$ and $\sqrt{3}$. See `ToyExamples.md` for the full
catalog of equations, slices, and known solutions.

The test loop at the end of `HybridDecomposition.m2` first checks this fixture
data algebraically, then passes each witness set to `componentRepresentation`.
It does not call `regeneration`, `decompose`, or `hybridDecomposition`. Run all
toy examples from `src` with

```sh
M2 -q --stop -e 'load "HybridDecomposition.m2"; testToyExamples();'
```

## Equidimensional benchmark test

`testSys2449EquidimensionalComponents()` exercises a larger system from
`ExampleSystems.m2`. The top-dimensional part of `Sys2449` has dimension
$d=8$, degree $4$, and exactly two irreducible components, of degrees $1$ and
$3$. The test generates eight deterministic random hyperplanes over
$\mathbb{Q}$, verifies that the resulting slice has dimension zero and degree
$4$, and uses NAG monodromy to group its four points into the two irreducible
components.

For each group, the test rebuilds a witness set with the rational slice, runs
`componentRepresentation`, solves the reconstructed $f(T)$, evaluates

$$
x_i=\frac{h_i(T)}{d_i f_{\mathrm{red}}'(T)},
$$

and compares the RUR-derived points bijectively with the sliced component at
tolerance $10^{-20}$. It also substitutes every recovered point into the
original `Sys2449` equations and all eight hyperplanes.

Run the toy and benchmark suites together with

```sh
M2 -q --stop -e 'load "HybridDecomposition.m2"; testExamples();'
```

`ExampleSystems.m2` remains available separately and contains exact `QQ`
versions of all nine polynomial systems from `assets/examples`, plus the
older circle-and-points fixture.  Each system is exported individually and
through `exampleSystems`, indexed by the names in `exampleSystemNames`.

Each exported entry contains `"Ring"`, `"Polynomials"`, `"Name"`, and the
original finite-field `"SourceCharacteristic"`. For example:

```m2
load "ExampleSystems.m2";
E = getExampleSystem "pseudosing_deg_2_nvars_6";
I = ideal(E#"Polynomials");
```

Several of the larger systems have high path counts, so decomposing them can
take a long time; the direct toy-example loop avoids that work.
