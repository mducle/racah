function E = racah_FtoE(F)
% Converts the slater integrals F0,F2,F4,F6 to the Racah parameters E0,E1,E2,E3

if isnumeric(F) & length(F)==4
  E = [( F(1) - 10*F(2) -  33*F(3) -  286*F(4) );   ...
       (        70*F(2) + 231*F(3) + 2002*F(4) )/9; ...
       (           F(2) -   3*F(3) +    7*F(4) )/9; ...
       (         5*F(2) +   6*F(3) -   91*F(4) )/3];
else
  error('Input must be a vector [F0 F2 F4 F6]');
end
