using Oscar
P, (A21, A31, A32, A41, A42, A43, B1, B2, B3, B4, C2, C3, C4) = PolynomialRing(GF(1484137), ["A21", "A31", "A32", "A41", "A42", "A43", "B1", "B2", "B3", "B4", "C2", "C3", "C4"]);
F = [B1+B2+B3+B4-1, 2*B2*C2+2*B3*C3+2*B4*C4-1, 3*B2*C2^2+3*B3*C3^2+3*B4*C4^2-1, 6*A32*B3*C2+6*A42*B4*C2+6*A43*B4*C3-1, 4*B2*C2^3+4*B3*C3^3+4*B4*C4^3-1, 8*A32*B3*C2*C3+8*A42*B4*C2*C4+8*A43*B4*C3*C4-1, 12*A32*B3*C2^2+12*A42*B4*C2^2+12*A43*B4*C3^2-1, 24*A32*A43*B4*C2-1, -A21+C2, -A31-A32+C3, -A41-A42-A43+C4];
sort!(F, by = p -> total_degree(p));
result_file = open("./results_rec_decomp/Sys3011_result.jl", "a+")
