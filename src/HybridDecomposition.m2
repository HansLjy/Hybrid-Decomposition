needsPackage "NumericalAlgebraicGeometry";
needsPackage "LLLBases";
load "ExampleSystems.m2";
load "NodalCubicToy.m2";

-----------------------------------------------------------------------------
-- Approximate integer kernels
-----------------------------------------------------------------------------

primitiveVector = c -> (
    nonzero := select(c, x -> x != 0);
    if #nonzero == 0 then return null;
    g := fold(nonzero / (x -> abs x), gcd);
    result := apply(c, x -> x // g);
    firstNonzero := first select(result, x -> x != 0);
    if firstNonzero < 0 then -result else result
    );

relativeRelationResidual = (M,c) -> max apply(numrows M, i -> (
    value := sum(numcols M, j -> M_(i,j) * c#j);
    scale := 1 + sum(numcols M, j -> abs(M_(i,j)) * abs(c#j));
    abs(value) / scale
    ));

relationsAtScale = (M,scale,tolerance,maxHeight) -> (
    m := numcols M;
    realRows := apply(numrows M, i ->
        apply(m, j -> round(scale * realPart M_(i,j))));
    imaginaryRows := apply(numrows M, i ->
        apply(m, j -> round(scale * imaginaryPart M_(i,j))));
    lattice := matrix(ZZ, entries id_(ZZ^m) | realRows | imaginaryRows);
    reduced := LLL lattice;
    candidates := {};
    for j from 0 to numcols reduced - 1 do (
        c := primitiveVector apply(m, i -> reduced_(i,j));
        if c =!= null
        and max(c / abs) <= maxHeight
        and relativeRelationResidual(M,c) <= tolerance
        then candidates = append(candidates,c);
        );
    unique candidates
    );

integerRelations = method(Options => {
    RelationScale => 10^60,
    RelationTolerance => 10^(-12),
    RelationMaxHeight => 10^30
    });

integerRelations Matrix := o -> M -> (
    if numrows M == 0 or numcols M == 0 then
        error "integerRelations expects a nonempty matrix";
    scale := o.RelationScale;
    if class scale =!= ZZ or scale <= 0 then
        error "RelationScale must be a positive integer";
    firstRelations := relationsAtScale(
        M,scale,o.RelationTolerance,o.RelationMaxHeight);
    secondScale := max(1000, floor(scale/100));
    secondRelations := relationsAtScale(
        M,secondScale,o.RelationTolerance,o.RelationMaxHeight);
    stable := select(firstRelations, c -> member(c,secondRelations));
    new HashTable from {
        "Relations" => stable,
        "Residuals" => apply(stable, c -> relativeRelationResidual(M,c)),
        "FirstScaleRelations" => firstRelations,
        "SecondScaleRelations" => secondRelations,
        "Scale" => scale,
        "SecondScale" => secondScale
        }
    );

relationOptions = (M,o) -> integerRelations(M,
    RelationScale => o.RelationScale,
    RelationTolerance => o.RelationTolerance,
    RelationMaxHeight => o.RelationMaxHeight);

-- Numerical witness data
-----------------------------------------------------------------------------

pointMultiplicity = p -> if p.cache.?Multiplicity then p.cache.Multiplicity else 1;

refineAndPreserveMultiplicities = (W,multiplicities,o) -> (
    oldPoints := points W;
    if #oldPoints != #multiplicities then error "multiplicity list has the wrong length";
    newPoints := if o.RefineWitnessPoints then apply(#oldPoints, i -> refine(oldPoints#i,
        Software => o.DecompositionSoftware,
        Bits => o.RefinementBits,
        ErrorTolerance => o.RefinementTolerance,
        ResidualTolerance => o.RefinementTolerance,
        Iterations => o.RefinementIterations))
    else oldPoints;
    scan(#newPoints, i -> newPoints#i.cache.Multiplicity = multiplicities#i);
    witnessSet(W#Equations,W#Slice,newPoints)
    );

prepareWitnessSet = (W,o) -> (
    multiplicities := toList(points W / pointMultiplicity);
    singularFlags := toList(points W / (p -> status p === Singular));
    refinedW := refineAndPreserveMultiplicities(W,multiplicities,o);
    (refinedW,multiplicities,singularFlags)
    );

linearFormCoefficients = (n,seed,bound,attempt) -> (
    c := apply(n, i ->
        ((seed + attempt + 19*(i+1) + 7*(i+1)^2) % (2*bound+1)) - bound);
    if all(c, x -> x == 0) then {1} | drop(c,1) else c
    );

sliceSolutionsWithMultiplicity = (W,singularFlags,o) -> (
    system := equations W | slice W;
    rawSolutions := solveSystem(system,
        Software => o.DecompositionSoftware,
        Precision => o.RefinementBits);
    clusteredSolutions := solutionsWithMultiplicity(rawSolutions,
        Tolerance => o.ClusterTolerance);
    multiplicities := toList apply(#(points W), i -> (
        p := (points W)#i;
        -- Total-degree homotopies may return duplicate paths even at a
        -- regular solution.  Such duplicates are not local multiplicity.
        if not singularFlags#i then return 1;
        matches := select(clusteredSolutions, q ->
            areEqual(p,q,Tolerance => o.ClusterTolerance));
        if #matches != 1 then
            error "could not match a witness point to the clustered slice solutions";
        pointMultiplicity first matches
        ));
    (multiplicities,rawSolutions,clusteredSolutions)
    );

chooseSeparatingLinearForm = (witnessPoints,multiplicities,o) -> (
    n := #coordinates first witnessPoints;
    answer := null;
    attempt := 0;
    while answer === null and attempt < o.LinearFormAttempts do (
        c := linearFormCoefficients(n,o.LinearFormSeed,
            o.LinearFormCoefficientBound,attempt);
        values := toList apply(witnessPoints, p ->
            sum(n, i -> c#i * (coordinates p)#i));
        expanded := flatten apply(#values, j ->
            apply(multiplicities#j, k -> point {{values#j}}));
        clusters := solutionsWithMultiplicity(expanded,
            Tolerance => o.ClusterTolerance);
        if #clusters == #witnessPoints then answer = (c,values,clusters);
        attempt = attempt + 1;
        );
    if answer === null then
        error "could not find an integer linear form separating the witness points";
    answer
    );

clusterMultiplicityAt = (beta,clusters,tolerance) -> (
    matches := select(clusters, p ->
        abs((coordinates p)#0-beta) <= tolerance*(1+abs beta));
    if #matches != 1 then error "could not match a projected root to its cluster";
    pointMultiplicity first matches
    );

-----------------------------------------------------------------------------
-- f and the coordinate numerators
-----------------------------------------------------------------------------

fallingFactorial = (r,k) -> if k == 0 then 1 else product(k, j -> r-j);

multiplicityMatrix = (values,multiplicities) -> (
    N := sum multiplicities;
    C := ring first values;
    rows := flatten apply(#values, j -> apply(multiplicities#j, k ->
        apply(N+1, r -> if r < k then 0_C
            else promote(fallingFactorial(r,k),C) * values#j^(r-k))));
    matrix rows
    );

selectFullDegreeRelation = (relations,N) -> (
    candidates := select(relations, c -> #c == N+1 and c#N != 0);
    if #candidates == 0 then null else fold(drop(candidates,1), first candidates,
        (c,best) -> if max(c / abs) < max(best / abs) then c else best)
    );

primitivePolynomial = f -> (
    T := (gens ring f)#0;
    d := first degree f;
    coeffs := apply(d+1, i -> coefficient(T^i,f));
    g := fold(select(coeffs, x -> x != 0) / (x -> abs x), gcd);
    q := sum(d+1, i -> (coeffs#i // g) * T^i);
    if leadCoefficient q < 0 then -q else q
    );

polynomialFromRelation = c -> (
    P := ZZ(monoid[Variables => 1, VariableBaseName => "T"]);
    T := (gens P)#0;
    primitivePolynomial sum(#c, i -> c#i*T^i)
    );

polynomialCoefficients = f -> (
    T := (gens ring f)#0;
    apply(first degree f + 1, i -> coefficient(T^i,f))
    );

evaluateDerivativeAt = (f,beta) -> (
    coeffs := polynomialCoefficients f;
    if #coeffs <= 1 then 0*beta
    else sum(toList(1..(#coeffs-1)), i -> i*coeffs#i*beta^(i-1))
    );

coordinateRelationMatrix = (values,witnessPoints,coordinateIndex,fSquareFree) -> (
    s := first degree fSquareFree;
    rows := apply(#values, j -> (
        beta := values#j;
        derivativeValue := evaluateDerivativeAt(fSquareFree,beta);
        apply(s, r -> beta^r) |
            {-((coordinates witnessPoints#j)#coordinateIndex) * derivativeValue}
        ));
    matrix rows
    );

coordinateRepresentation = (M,P,o) -> (
    s := numcols M - 1;
    data := relationOptions(M,o);
    candidates := select(data#"Relations", c -> #c == s+1 and c#s != 0);
    if #candidates == 0 then return new HashTable from {
        "Reconstructed" => false, "Denominator" => null,
        "Numerator" => null, "RelationData" => data
        };
    c := fold(drop(candidates,1), first candidates,
        (candidate,best) -> if max(candidate / abs) < max(best / abs)
            then candidate else best);
    if c#s < 0 then c = -c;
    denominator := c#s;
    T := (gens P)#0;
    numerator := sum(s, r -> c#r*T^r);
    new HashTable from {
        "Reconstructed" => true,
        "Denominator" => denominator,
        "Numerator" => numerator,
        "Relation" => c,
        "Residual" => relativeRelationResidual(M,c),
        "RelationData" => data
        }
    );

componentRepresentation = method(Options => {
    RelationScale => 10^60, RelationTolerance => 10^(-12),
    RelationMaxHeight => 10^30, RefineWitnessPoints => true,
    RefinementBits => 300, RefinementTolerance => 10^(-70),
    RefinementIterations => 60, LinearFormSeed => 29,
    LinearFormCoefficientBound => 13, LinearFormAttempts => 100,
    ClusterTolerance => 10^(-6), DecompositionSoftware => M2engine,
    PrintResults => true
    });

componentRepresentation (WitnessSet,PolynomialRing) := o -> (W,inputRing) -> (
    (refinedW,cachedMultiplicities,singularFlags) := prepareWitnessSet(W,o);
    witnessPoints := toList points refinedW;
    (sliceMultiplicities,rawSliceSolutions,clusteredSliceSolutions) :=
        sliceSolutionsWithMultiplicity(refinedW,singularFlags,o);
    (linearCoefficients,values,clusters) := chooseSeparatingLinearForm(
        witnessPoints,sliceMultiplicities,o);
    multiplicities := toList apply(values, beta ->
        clusterMultiplicityAt(beta,clusters,o.ClusterTolerance));
    N := sum multiplicities;
    Mf := multiplicityMatrix(values,multiplicities);
    fRelationData := relationOptions(Mf,o);
    fRelation := selectFullDegreeRelation(fRelationData#"Relations",N);
    exactLinearForm := sum(numgens inputRing, i ->
        linearCoefficients#i * (gens inputRing)#i);

    baseRecord := new HashTable from {
        "Dimension" => dim refinedW, "Degree" => degree refinedW,
        "Hyperplanes" => slice refinedW,
        "WitnessSet" => refinedW, "WitnessPoints" => witnessPoints,
        "Multiplicities" => multiplicities, "LinearForm" => exactLinearForm,
        "CachedMultiplicities" => cachedMultiplicities,
        "SingularWitnessPoints" => singularFlags,
        "RawSliceSolutions" => rawSliceSolutions,
        "ClusteredSliceSolutions" => clusteredSliceSolutions,
        "LinearFormCoefficients" => linearCoefficients,
        "ProjectedValues" => values, "FRelationData" => fRelationData
        };
    if fRelation === null then return merge(baseRecord,new HashTable from {
        "f" => null, "SquareFreePart" => null,
        "CoordinateRepresentations" => {}, "Reconstructed" => false,
        "Failure" => "no stable full-degree integer relation for f"
        },(a,b)->b);

    if fRelation#N < 0 then fRelation = -fRelation;
    f := polynomialFromRelation fRelation;
    T := (gens ring f)#0;
    fSquareFree := primitivePolynomial(f // gcd(f,diff(T,f)));
    if first degree fSquareFree != #witnessPoints then return merge(baseRecord,
        new HashTable from {
            "f" => f, "SquareFreePart" => fSquareFree,
            "CoordinateRepresentations" => {}, "Reconstructed" => false,
            "Failure" => "squarefree degree differs from the number of distinct witness points"
            },(a,b)->b);

    coordinateData := apply(numgens inputRing, i -> coordinateRepresentation(
        coordinateRelationMatrix(values,witnessPoints,i,fSquareFree),ring f,o));
    success := all(coordinateData, h -> h#"Reconstructed");
    if o.PrintResults then (
        print("dimension " | toString(dim refinedW) | ", degree " | toString(degree refinedW));
        print("f = " | toString f);
        if success then scan(#coordinateData, i -> print(
            "g_" | toString(i+1) | " = (" |
            toString(coordinateData#i#"Numerator") | ")/" |
            toString(coordinateData#i#"Denominator")))
        else print "one or more coordinate relations could not be reconstructed";
        );
    merge(baseRecord,new HashTable from {
        "f" => f, "SquareFreePart" => fSquareFree,
        "CoordinateRepresentations" => coordinateData,
        "Denominators" => coordinateData / (h -> h#"Denominator"),
        "NumeratorPolynomials" => coordinateData / (h -> h#"Numerator"),
        "Reconstructed" => success,
        "Failure" => if success then null else "a coordinate relation failed",
        "FRelation" => fRelation,
        "FResidual" => relativeRelationResidual(Mf,fRelation)
        },(a,b)->b)
    );

-----------------------------------------------------------------------------
-- Full numerical decomposition
-----------------------------------------------------------------------------

hybridDecomposition = method(Options => {
    RelationScale => 10^60, RelationTolerance => 10^(-12),
    RelationMaxHeight => 10^30, RefineWitnessPoints => true,
    RefinementBits => 300, RefinementTolerance => 10^(-70),
    RefinementIterations => 60, LinearFormSeed => 29,
    LinearFormCoefficientBound => 13, LinearFormAttempts => 100,
    ClusterTolerance => 10^(-6), SeparateComponents => true,
    DecompositionSoftware => M2engine, DecompositionOutput => Singular,
    PrintResults => true
    });

hybridDecomposition Ideal := o -> I -> (
    R := ring I;
    if coefficientRing R =!= QQ then
        error "expected an ideal in a polynomial ring over QQ";
    V := regeneration(I_*,Software => o.DecompositionSoftware,
        Output => o.DecompositionOutput);
    dimensions := sort select(keys V, d -> class d === ZZ);
    byDimension := new MutableHashTable;
    scan(dimensions, d -> (
        Ws := if o.SeparateComponents then
            flatten apply(V#d, W -> decompose W) else V#d;
        byDimension#d = apply(#Ws, componentIndex -> (
            record := componentRepresentation(Ws#componentIndex,R,
                RelationScale => o.RelationScale,
                RelationTolerance => o.RelationTolerance,
                RelationMaxHeight => o.RelationMaxHeight,
                RefineWitnessPoints => o.RefineWitnessPoints,
                RefinementBits => o.RefinementBits,
                RefinementTolerance => o.RefinementTolerance,
                RefinementIterations => o.RefinementIterations,
                LinearFormSeed => o.LinearFormSeed + d + componentIndex,
                LinearFormCoefficientBound => o.LinearFormCoefficientBound,
                LinearFormAttempts => o.LinearFormAttempts,
                ClusterTolerance => o.ClusterTolerance,
                DecompositionSoftware => o.DecompositionSoftware,
                PrintResults => o.PrintResults);
            merge(record,new HashTable from {"ComponentIndex"=>componentIndex},(a,b)->b)
            ));
        ));
    new HashTable from {
        "InputIdeal" => I, "NumericalVariety" => V,
        "Dimensions" => dimensions, "ByDimension" => byDimension
        }
    );

-- Direct witness-set test on the rational nodal-cubic toy
-----------------------------------------------------------------------------

testNodalCubicToyWitnessSet = () -> (
    curveEquation := NodalCubicToyCurveIdeal_0;
    nodeSubstitution := {nodalX=>0_QQ,nodalY=>0_QQ};
    assert(sub(curveEquation,nodeSubstitution) == 0);
    assert(sub(diff(nodalX,curveEquation),nodeSubstitution) == 0);
    assert(sub(diff(nodalY,curveEquation),nodeSubstitution) == 0);
    assert(all(NodalCubicToyIsolatedCoordinates, p ->
        sub(curveEquation,{nodalX=>p#0,nodalY=>p#1}) != 0));
    assert(all(NodalCubicToyIsolatedCoordinates, p ->
        all(NodalCubicToyIdeal_*, f ->
            sub(f,{nodalX=>p#0,nodalY=>p#1}) == 0)));
    assert(all(NodalCubicToyIsolatedCoordinates, p ->
        sub(NodalCubicToySliceIdeal_0,
            {nodalX=>p#0,nodalY=>p#1}) != 0));
    assert(dim NodalCubicToySliceIntersectionIdeal == 0);
    assert(degree NodalCubicToySliceIntersectionIdeal == 3);

    assert(dim NodalCubicToyWitnessSet == 1);
    assert(degree NodalCubicToyWitnessSet == 3);
    assert(all(NodalCubicToyWitnessCoordinates, p ->
        sub(curveEquation,{nodalX=>p#0,nodalY=>p#1}) == 0));
    assert(all(NodalCubicToyWitnessCoordinates, p ->
        all(NodalCubicToyIdeal_*, f ->
            sub(f,{nodalX=>p#0,nodalY=>p#1}) == 0)));
    assert(all(NodalCubicToyWitnessCoordinates, p ->
        sub(NodalCubicToySliceIdeal_0,{nodalX=>p#0,nodalY=>p#1}) == 0));

    expectedR := product(NodalCubicToyWitnessCoordinates, p ->
        nodalU0-p#0*nodalU1-p#1*nodalU2);
    assert(NodalCubicToyR == expectedR);

    setRandomSeed 0;
    print "testing the known nodal-cubic witness set";
    result := componentRepresentation(
        NodalCubicToyWitnessSet,
        NodalCubicToyRing,
        RefineWitnessPoints => false,
        PrintResults => true);
    assert(result#"Reconstructed");
    assert(result#"Dimension" == 1);
    assert(result#"Degree" == 3);
    assert(result#"WitnessSet" === NodalCubicToyWitnessSet);
    assert(#(result#"ClusteredSliceSolutions") == 3);

    expectedF := nodalCubicToySpecializeR(
        result#"LinearFormCoefficients");
    expectedRing := ring expectedF;
    reconstructedF := (map(
        expectedRing,
        ring result#"f",
        {(gens expectedRing)#0}))(result#"f");
    assert(reconstructedF / leadCoefficient reconstructedF ==
        expectedF / leadCoefficient expectedF);
    print "nodal-cubic direct witness-set test completed";
    result
    );

end;
