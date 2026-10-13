needsPackage "NumericalAlgebraicGeometry";

-----------------------------------------------------------------------------
-- A rational witness set for a nodal cubic together with isolated points
-----------------------------------------------------------------------------

-- The curve y^2 = x^2(x+1) has a node at (0,0).  The full toy variety is
-- this curve union three isolated rational points away from the curve.
NodalCubicToyRing = QQ[nodalX,nodalY];
NodalCubicToyCurveIdeal = ideal(
    nodalY^2 - nodalX^2*(nodalX+1));
NodalCubicToyIsolatedCoordinates = {
    {0_QQ,2_QQ},
    {1_QQ,2_QQ},
    {2_QQ,2_QQ}
    };
NodalCubicToyIsolatedIdeals = apply(
    NodalCubicToyIsolatedCoordinates, p ->
        ideal(nodalX-p#0,nodalY-p#1));
NodalCubicToyIdeal = fold(
    NodalCubicToyIsolatedIdeals,
    NodalCubicToyCurveIdeal,
    (pointIdeal,unionIdeal) -> intersect(pointIdeal,unionIdeal));

-- The rational line 5y-18x+24=0 meets the cubic in the three distinct
-- rational points below.  They arise from the standard parametrization
-- (x,y)=(s^2-1,s(s^2-1)) at s=2,3,-7/5.
NodalCubicToySliceIdeal = ideal(5*nodalY-18*nodalX+24);
NodalCubicToySliceIntersectionIdeal =
    NodalCubicToyIdeal + NodalCubicToySliceIdeal;
NodalCubicToyWitnessCoordinates = {
    {3_QQ,6_QQ},
    {8_QQ,24_QQ},
    {24/25,-168/125}
    };

-- Store the known exact coordinates at ample numerical precision.  No
-- numerical decomposition or path tracking is needed to construct this
-- witness set.
NodalCubicToyCC = CC_300;
NodalCubicToyWitnessPoints = apply(
    NodalCubicToyWitnessCoordinates, p ->
        point {apply(p, q -> promote(q,NodalCubicToyCC))});
-- This is the dimension-one component's witness set.  The three isolated
-- components of NodalCubicToyIdeal do not lie on the chosen slice.
NodalCubicToyWitnessSet = witnessSet(
    NodalCubicToyCurveIdeal,
    NodalCubicToySliceIdeal,
    NodalCubicToyWitnessPoints);
declareIrreducible NodalCubicToyWitnessSet;

-- Closed-form Chow polynomial from the progress report:
-- R(U0,U1,U2) = product_alpha (U0-alpha_1 U1-alpha_2 U2).
NodalCubicToyChowRing = QQ[nodalU0,nodalU1,nodalU2];
NodalCubicToyR =
    (nodalU0-3*nodalU1-6*nodalU2) *
    (nodalU0-8*nodalU1-24*nodalU2) *
    (nodalU0-(24/25)*nodalU1+(168/125)*nodalU2);

nodalCubicToySpecializeR = linearCoefficients -> (
    if #linearCoefficients != 2 then
        error "the nodal-cubic toy linear form needs two coefficients";
    P := QQ[expectedT];
    phi := map(P,NodalCubicToyChowRing,{
        expectedT,
        promote(linearCoefficients#0,P),
        promote(linearCoefficients#1,P)
        });
    phi NodalCubicToyR
    );
