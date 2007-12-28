function out = racahW(A)
% racahW - Calculates the value of the Racah W function: W(abcd;ef)
% 
% Inputs:  A is a 6-component vector [a b c d e f]
% Outputs: W = W(abcd;ef) is a double scalar.
%
% The Racah W function is given in Racah, Phys. Rev., vol 62, pp438 (1942) eqn (36).

% By Duc Le - duc.le@ucl.ac.uk 2006

% The formula for the Racah W function is:
%
%                                        1/2 
% W(abcd;ef) = [T(abe)T(cde)T(acf)T(bdf)]     
% 
%               ---                     (a + b + c + d + 1 - z)!
%             * >    ------------------------------------------------------------------
%               ---z z!(a+b-e-z)!(c+d-e-z)!(a+c-f-z)!(b+d-f-z)!(e+f-a-d+z)!(e+f-b-c+z)!
%
% where the triangle functions are:
%
%          (a+b-c)!(a-b+c)!(-a+b+c)!
% T(abc) = -------------------------
%              (a + b + c + 1)!
%
% The W functions are only defined for integer or half integer arguments that satisfy
% the selection rule that the four triads (abe), (cde), (acf), and (bdf) has an integral
% sum and satisfy the triangular inequality |a-b| <= c <= |a+b|. Otherwise W is zero.

[szrow, szcol] = size(A);
if (szrow ~= 1) | (szcol ~= 6)
  error('Input requires a 1x6 matrix')
end

a = A(1); b = A(2); c = A(3);
d = A(4); e = A(5); f = A(6);

% Defines a 4x3 matrix to store the elements of the triangular coefficients
m = [a b e;
     c d e;
     a c f;
     b d f];

% Selection rules:
% The triads (a b e), (c d e), (a c f), (b d f)
% 1. All satisfy the triangular inequality: |a-b| <= c <= a+b  
if ~isequal(zeros(4,1), (m(:,3) < abs(m(:,1)-m(:,2))) | (m(:,3) > m(:,1)+m(:,2)))
  out = 0; return;
% NB relational operators compares element by element and returns a matrix of 1 or 0 
% 2. Elements of each triad sum to an integer
%elseif ~isequal(zeros(1,4), double(int8(sum(m'))) ./ sum(m') - 1)
elseif ~isequal(zeros(1,4), mod(sum(m'),1))
  out = 0; return;
end

% Calculates the triangular coefficients
tr = prod(factorial([m(:,1)+m(:,2)-m(:,3) m(:,1)-m(:,2)+m(:,3) -m(:,1)+m(:,2)+m(:,3)])');
tb = factorial(m(:,1)+m(:,2)+m(:,3)+1);
tr = tr(:); tb = tb(:);

% Defines a vector to store individual factorial arguments of f(t)
F = [a+b-e c+d-e a+c-f b+d-f e+f-a-d e+f-b-c];

% Calculates the sum over z
%sum_z = 0;
%for z = max([0 -F(5:6)]):min([F(1:4) a+b+c+d+1])
%    sum_z = sum_z + (-1)^z * factorial(a+b+c+d+1-z) / prod(factorial([z F(1:4)-z F(5:6)+z]));
%end
z = max([0 -F(5:6)]):min([F(1:4) a+b+c+d+1]);
out = sqrt(prod(tr./tb)) * ...% sum_z;
      sum( (-1).^z .* factorial(a+b+c+d+1-z) ./ prod(factorial([z;F(1)-z;F(2)-z;F(3)-z;F(4)-z;F(5)+z;F(6)+z])) );
