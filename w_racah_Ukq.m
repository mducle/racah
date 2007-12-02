function [Ukq,states] = racah_Ukq(n,l,k,indq,statesLS,Jind)
% Calculates the matrix Ukq, after Elliot, Judd, and Runciman.

icfact = [-sqrt(15/7)/2       sqrt(11/14)        -sqrt(429/7)/10  ];

if ~exist('indq')
  indq = 1:(2*k+1);
else
  indq = indq + k+1;
end

if ~exist('statesLS')
  statesLS = racah_states(n,l);
end

if ~exist('Jind'); Jind_flag = 1; end

redmat = racah_Umat(n,l,k,statesLS);

index = 0;
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2});
  S = statesLS{i}{1};
  Jmin = abs(L-S); Jmax = L+S;
  if exist('Jind_flag'); Jind = Jmin:Jmax; end
  for iJ = Jind
    for mJ = -iJ:iJ
      index = index + 1;
      states{index} = {S L statesLS{i}{3:4} iJ mJ i};
    end
  end
end

num_states = length(states);

if size(indq~=1)
  for iq = indq; Ukq{iq} = zeros(num_states); end;
else
  Ukq{iq} = zeros(num_states);
end
for i = 1:num_states
  L = states{i}{2}; S = states{i}{1}; J = states{i}{5}; Jz = states{i}{6}; irm = states{i}{7};
  for j = 1:num_states 
    Lp = states{j}{2}; Sp = states{j}{1}; Jp = states{j}{5}; Jzp = states{j}{6}; irmp = states{j}{7};
    rm = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]) * redmat(irm,irmp);
    for iq = indq
      q = iq-1-k;
      % NB. wigner symbol (j1j2jm|j1j2m1m2) = (-1)^(j1-j2-m) * sqrt(2*j+1) * [j1 j2 j; m1 m2 m]
      if S==Sp
        Ukq{iq}(i,j) = (-1)^(J+Jz+k+q) * (1/sqrt(2*k+1)) * rm * wigner([J Jp -Jz Jzp J Jp k -q]) / icfact(k/2);
      end
    end
  end
end

if length(indq)==1
  Ukq = Ukq{indq};
end
