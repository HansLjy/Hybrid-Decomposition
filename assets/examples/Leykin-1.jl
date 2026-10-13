using Oscar
P, (x4, x3, x2, x1, a4, a3, a2, a1) = PolynomialRing(GF(1909907), ["x4", "x3", "x2", "x1", "a4", "a3", "a2", "a1"]);
F = [(x1^2*x3-x1*x2^2)*x4-x1*x2*x3^2+x2^3*x3, 5*a1*x1+a2*x2-3*a3*x3-7*a4*x4, -4*a4*x2*x4-5*a4*x3^2+(-5*a2*x1-6*a3*x2)*x3+2*a2*x2^2, -2*a4*x1*x4+(2*a3*x1+5*a4*x2)*x3+a2*x1*x2, a2*x1^2+2*a3*x1*x2+a4*x1*x3+2*a4*x2^2, (5*a2*x1+10*a3*x2+7*a4*x3)*x4-2*a3*x3^2+4*a2*x2*x3+10*a1*x2*x1];
sort!(F, by = p -> total_degree(p));
result_file = open("./results_rec_decomp/Leykin-1_result.jl", "a+")
