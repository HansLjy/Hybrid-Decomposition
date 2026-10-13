using Oscar
P, (x1, x2, x3, y1, y2, y3)  = PolynomialRing(GF(1063441), ["x1", "x2", "x3", "y1", "y2", "y3"]);
F = [-7*x1^2+22*x1*x2-55*x1*x3+87*x2^2-56*x2*x3-62*x3^2-94*x1+97*x3-73, -4*x1^2-83*x1*x2-10*x1*x3-82*x2^2+80*x2*x3+71*x3^2+62*x1-44*x2-17*x3-75, -7*y1^2+22*y1*y2-55*y1*y3+87*y2^2-56*y2*y3-62*y3^2-94*y1+97*y3-73, -4*y1^2-83*y1*y2-10*y1*y3-82*y2^2+80*y2*y3+71*y3^2+62*y1-44*y2-17*y3-75, x3-y3, x2-y2];
sort!(F, by = p -> total_degree(p));
result_file = open("./results_rec_decomp/pseudosing_deg_2_nvars_6_result.jl", "a+")
