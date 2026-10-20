needsPackage "NumericalAlgebraicGeometry";
needsPackage "LLLBases";
load "ExampleSystems.m2";
load "ToyExamples.m2";

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
        print("linear form = " | toString exactLinearForm);
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

-- Direct witness-set tests on the explicitly solved toy examples
-----------------------------------------------------------------------------

toyPolynomialResidual = (f,p) -> (
    R := ring f;
    z := coordinates p;
    abs sub(f,apply(numgens R, i -> (gens R)#i => z#i))
    );

evaluateUnivariateAt = (f,beta) -> (
    coefficients := polynomialCoefficients f;
    sum(#coefficients, i -> coefficients#i*beta^i)
    );

solveRUR = result -> (
    if not result#"Reconstructed" then
        error "cannot solve an RUR that was not reconstructed";
    f := result#"f";
    fRed := result#"SquareFreePart";
    coordinateData := result#"CoordinateRepresentations";
    projectedRoots := roots(f,Precision => 300);
    apply(projectedRoots, beta -> (
        derivativeValue := evaluateDerivativeAt(fRed,beta);
        if abs derivativeValue < 10^(-70) then
            error "encountered a multiple root of the squarefree part";
        point {apply(coordinateData, coordinateRecord ->
            evaluateUnivariateAt(coordinateRecord#"Numerator",beta) /
                (coordinateRecord#"Denominator"*derivativeValue))}
        ))
    );

samePointSets = (leftPoints,rightPoints,tolerance) -> (
    #leftPoints == #rightPoints
    and all(leftPoints, p -> #select(rightPoints, q ->
        areEqual(p,q,Tolerance => tolerance)) == 1)
    and all(rightPoints, q -> #select(leftPoints, p ->
        areEqual(p,q,Tolerance => tolerance)) == 1)
    );

printToySolutions = (exactDescriptions,numericalPoints) -> (
    oldPrintingPrecision := printingPrecision;
    printingPrecision = 20;
    print("exact solutions = " | toString exactDescriptions);
    print("numerical RUR solutions = " |
        toString apply(numericalPoints, coordinates));
    printingPrecision = oldPrintingPrecision;
    );

testToyExample = name -> (
    example := toyExamples#name;
    n := example#"n";
    d := example#"d";
    R := example#"Ring";
    I := example#"Ideal";
    H := example#"Hyperplanes";
    knownSolutions := example#"KnownSolutions";

    print("testing toy example " | name |
        " with (n,d) = (" | toString n | "," | toString d | ")");
    assert(coefficientRing R === QQ);
    assert(numgens R == n);
    assert(dim I == d);
    assert(numgens H == d);
    assert(all(H_*, h -> first degree h == 1));
    assert(dim(I+H) == 0);
    assert(degree(I+H) == #knownSolutions);
    assert(all(knownSolutions, p ->
        all(I_* | H_*, f -> toyPolynomialResidual(f,p) < 10^(-70))));

    result := componentRepresentation(
        example#"WitnessSet",
        R,
        RefineWitnessPoints => false,
        PrintResults => true);
    assert(result#"Reconstructed");
    assert(result#"Dimension" == d);
    assert(result#"Degree" == #knownSolutions);
    assert(#(result#"WitnessPoints") == #knownSolutions);
    rurSolutions := solveRUR result;
    printToySolutions(example#"KnownSolutionDescriptions",rurSolutions);
    assert(samePointSets(rurSolutions,knownSolutions,10^(-20)));
    print("toy example " | name | " completed");
    result
    );

testToyExamples = () -> (
    setRandomSeed 0;
    apply(toyExampleNames, testToyExample)
    );

-- Keep the old focused entry point for callers that only want the nodal cubic.
testNodalCubicToyWitnessSet = () -> (
    setRandomSeed 0;
    testToyExample "nodal-cubic-rational-slice"
    );

-----------------------------------------------------------------------------
-- Equidimensional-component test from ExampleSystems.m2
-----------------------------------------------------------------------------

randomRationalHyperplanes = (R,d,seed,coefficientBound) -> (
    setRandomSeed seed;
    variables := gens R;
    ideal apply(d, j ->
        sum(numgens R, i ->
            random(-coefficientBound,coefficientBound)*variables#i)
        + random(-coefficientBound,coefficientBound))
    );

squareComponentIdeal = J -> (
    K := ideal take(J_*,codim J);
    if not isSubset(J,K) or not isSubset(K,J) then
        error "the selected square component equations do not define the component";
    K
    );

solveAndRefineSquareSlice = (K,H) -> (
    F := K_* | H_*;
    rawSolutions := solveSystem(F,Software => M2engine,Precision => 300);
    system := polySystem F;
    apply(rawSolutions, p -> refine(system,p,
        Software => M2engine,
        Bits => 300,
        ErrorTolerance => 10^(-70),
        ResidualTolerance => 10^(-70),
        Iterations => 60))
    );

complexWitnessSetForDecomposition = (I,H,slicePoints) -> (
    R := ring I;
    C := CC_100(monoid[
        Variables => numgens R,
        VariableBaseName => "decompositionVariable"]);
    phi := map(C,R,gens C);
    complexI := ideal apply(I_*, f -> phi f);
    complexH := ideal apply(H_*, f -> phi f);
    witnessSet(complexI,complexH,slicePoints)
    );

printNumericalPointSet = (label,pointSet) -> (
    oldPrintingPrecision := printingPrecision;
    printingPrecision = 20;
    print(label | " = " | toString apply(pointSet, coordinates));
    printingPrecision = oldPrintingPrecision;
    );

testSys2449EquidimensionalComponents = () -> (
    example := getExampleSystem "Sys2449";
    R := example#"Ring";
    I := ideal(example#"Polynomials");
    d := dim I;
    topComponents := select(minimalPrimes I, J -> dim J == d);

    assert(numgens R == 11);
    assert(d == 8);
    assert(#topComponents == 2);
    assert(sort apply(topComponents, degree) == {1,3});
    assert(sum(topComponents / degree) == degree I);

    H := randomRationalHyperplanes(R,d,2449,5);
    assert(coefficientRing ring H === QQ);
    assert(numgens H == d);
    assert(all(H_*, h -> first degree h == 1));
    assert(dim(I+H) == 0);
    assert(degree(I+H) == 4);

    print "testing Sys2449 top-dimensional equidimensional part";
    print("dimension = " | toString d |
        ", irreducible component degrees = " |
        toString sort apply(topComponents, degree));
    print("random rational hyperplanes = " | toString H_*);

    initialPointsByComponent := apply(topComponents, J -> (
        K := squareComponentIdeal J;
        slicePoints := solveAndRefineSquareSlice(K,H);
        assert(#slicePoints == degree J);
        assert(all(slicePoints, p -> all(J_* | H_*, f ->
            toyPolynomialResidual(f,p) < 10^(-70))));
        slicePoints
        ));
    initialSlicePoints := flatten initialPointsByComponent;
    assert(#initialSlicePoints == 4);
    assert(all(initialSlicePoints, p -> all(I_* | H_*, f ->
        toyPolynomialResidual(f,p) < 10^(-70))));

    -- NAG's monodromy decomposition needs a complex coefficient ring for
    -- its temporary slices.  The hyperplanes themselves remain rational,
    -- and rational witness sets are rebuilt below for reconstruction.
    complexWitnessSet := complexWitnessSetForDecomposition(
        I,H,initialSlicePoints);
    setRandomSeed 3011;
    numericalComponents := decompose complexWitnessSet;
    assert(#numericalComponents == 2);
    assert(sort apply(numericalComponents, degree) == {1,3});
    print("numerically grouped component degrees = " |
        toString sort apply(numericalComponents, degree));

    componentResults := apply(#numericalComponents, componentIndex -> (
        numericalComponent := numericalComponents#componentIndex;
        componentDegree := degree numericalComponent;
        matchingComponents := select(topComponents, J ->
            degree J == componentDegree);
        assert(#matchingComponents == 1);
        J := first matchingComponents;
        K := squareComponentIdeal J;
        system := polySystem(K_* | H_*);
        refinedPoints := apply(toList points numericalComponent, p ->
            refine(system,p,
                Software => M2engine,
                Bits => 300,
                ErrorTolerance => 10^(-70),
                ResidualTolerance => 10^(-70),
                Iterations => 60));
        assert(#refinedPoints == componentDegree);
        assert(all(refinedPoints, p -> all(I_* | H_*, f ->
            toyPolynomialResidual(f,p) < 10^(-70))));

        rationalWitnessSet := witnessSet(K,H,refinedPoints);
        declareIrreducible rationalWitnessSet;
        print("reconstructing Sys2449 component " |
            toString(componentIndex+1) | " of degree " |
            toString componentDegree);
        result := componentRepresentation(
            rationalWitnessSet,
            R,
            RefineWitnessPoints => false,
            LinearFormSeed => 29 + componentIndex,
            PrintResults => true);
        assert(result#"Reconstructed");
        assert(result#"Dimension" == d);
        assert(result#"Degree" == componentDegree);

        rurSolutions := solveRUR result;
        printNumericalPointSet("numerical slice solutions",refinedPoints);
        printNumericalPointSet("numerical RUR solutions",rurSolutions);
        assert(samePointSets(rurSolutions,refinedPoints,10^(-20)));
        assert(all(rurSolutions, p -> all(I_* | H_*, f ->
            toyPolynomialResidual(f,p) < 10^(-20))));
        result
        ));
    print "Sys2449 equidimensional-component test completed";
    new HashTable from {
        "Example" => example,
        "Ideal" => I,
        "Dimension" => d,
        "Hyperplanes" => H,
        "NumericalComponents" => numericalComponents,
        "ComponentResults" => componentResults
        }
    );

testExamples = () -> new HashTable from {
    "ToyExamples" => testToyExamples(),
    "Sys2449" => testSys2449EquidimensionalComponents()
    };

end;
