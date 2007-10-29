function g = racah_g(U)
% Calculates the Racah function g(U) or g(W), the eigenvalues of G(G2) or G(R7), after Racah IV

if isvector(U) & length(U)==2
  g = ( U(1)^2 + U(1)*U(2) + U(2)^2 + 5*U(1) + 4*U(2) ) / 12;
elseif isvector(U) & length(U)==3
  W = U;
  g = ( W(1)*(W(1)+5) + W(2)*(W(2)+3) + W(3)*(W(3)+1) ) / 10; 
else
  error('U must be a vector of length 2');
end
