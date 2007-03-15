function g = racah_g(U)
% Calculates the Racah function g(U), the eigenvalues of G(G2), after Racah IV

if isvector(U) & length(U)==2
  g = (U(1)^2 + U(1)*U(2) + U(2)^2 + 5*U(1) + 4*U(2))/12;
else
  error('U must be a vector of length 2');
end
