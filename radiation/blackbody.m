function I = blackbody(T,nu)
C1 = 3.741771790075259e-16;
C2 = 0.014387741858429;
I = 1.0 / pi * C1 * (nu*100).^3 ./ (exp(C2*nu*100./T)-1);
end

