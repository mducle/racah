function Hcfic = cfic_hmltn(n,B2,B4,B6,LSstate,J)
% cfic_hmltn - calculates the CF hamiltonian using Racah operators for a specific |vULSJ> state
%
% Syntax:  Hcfic = cfic_hmltn(n,B2,B4,B6,LSstate,J)
%          Hcfic = cfic_hmltn(n,B,LSstate)
%
% Inputs:  n       - scalar - number of f electrons in the unfilled valence shell
%          B       - cell   - {B2 B4 B6} cell array of CF parameters
%          B2      - vector - a length 5 vector of the second order CF parameters
%          B4      - vector - a length 9 vector of the 4th order CF parameters
%          B6      - vector - a length 13 vector of the 6th order CF parameters
%          LSstate - cell   - either {S L v [U1 U2]} or {S L v [U1 U2] J}
%          J       - scalar - total angular momentum quantum number.
%
% Outputs: Hcfic   - matrix - a (2J+1)x(2J+1) matrix representing the CF operator
%
% NB. If LSstate is a length 5 cell, you must also specify J. 

icfact = [-0.73192505471140   0.88640526042792   -0.78285192907544];
icfact = [1 1 1];
icfact = [-sqrt(15/28)        sqrt(11/14)        -sqrt(429/700)   ];

% Parses inputs

if iscell(B4) 
  if length(B4)==4
    if isscalar(B6)
      LSstate = B4;
      J = B6;
    else
      error('If LSstate is a length 4 cell, you must also specify J');
    end
  elseif length(B4)==5
    J = LSstate{5};
    LSstate = LSstate(1:4);
  else
    error('LSstate is not a length 4 or length 5 cell');
  end
end

if iscell(B2) && length(B2)==3
  if length(B2{2})==9;  B4 = B2{2}; else; error('B4 must be length 9');  end
  if length(B2{3})==13; B6 = B2{3}; else; error('B6 must be length 13'); end
  if length(B2{1})==5;  B2 = B2{1}; else; error('B2 must be length 5');  end
elseif isvector(B2)
  if length(B2)~=5; error('B2 must be length 5'); end
else
  error('B must be a length 3 cell array');
end

if iscell(LSstate) 
  if length(LSstate)==4 
    if ~exist('J') || ~isscalar(J)
      error('If LSstate is a length 4 cell, you must also specify a scalar J');
    end
  elseif length(LSstate)==5
    J = LSstate{5};
    LSstate = LSstate(1:4);
  else
    error('LSstate must be a length 4 or 5 cell array');
  end
else
  error('LSstate must be a cell array');
end

% Determines the pairs of k,q which are non-zero
iq2 = find(B2);
iq4 = find(B4);
iq6 = find(B6);

indq = {iq2 iq4 iq6};
indk = find([sum(iq2) sum(iq4) sum(iq6)]); 

B = {B2 B4 B6};

L = LSstate{2}; S = LSstate{1};
if isstr(L); L = racah_lconv(L); end
Lp = L; Sp = S; Jp = J;

Hcfic = zeros(2*J+1);

for ik = indk
  k = 2*ik;
  % Finds out the reduced matrix elements <vULS|U^k|vULS>
  redmat = racah_Umat(n,3,k,{LSstate});
  rm = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]) * redmat;

  for i = 1:(2*J+1)
    Jz = i-1-J;
    for j = 1:(2*J+1)
      Jzp = j-1-J;
      
      for iq = indq{ik}
        q = iq-1-k;
        Hcfic(i,j) = Hcfic(i,j) + (-1)^(J+Jz+k+q) * (1/sqrt(2*k+1)) * rm * wigner([J Jp -Jz Jzp J Jp k -q]) * B{ik}(iq) / icfact(ik);
      end  % q loop

    end    % j loop
  end      % i loop

end        % k loop
