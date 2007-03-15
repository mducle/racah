function [Ukq,states] = racah_Ukq(n,l,k,indq)
% Calculates the matrix Ukq, after Elliot, Judd, and Runciman.

if mod(k,2)~=0 || k>6 || k<2
  error('Rank k must be 2,4,6');
end

if ~exist('indq')
  indq = 1:(2*k+1);
else
  indq = indq + k+1;
end

%redmat = racah_Umat(n,l,k);
load(sprintf('redmat%1g.mat',k));

statesLS = racah_states(n,l);

index = 0;
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2});
  S = statesLS{i}{1};
  Jmin = abs(L-S); Jmax = L+S;
  %iredmat = 0;
  for iJ = Jmin:Jmax
    %iredmat = iredmat + 1;
    for mJ = -iJ:iJ
      index = index + 1;
      %states{index} = {S L statesLS{i}{3:4} iJ mJ iredmat};
      states{index} = {S L statesLS{i}{3:4} iJ mJ i};
    end
%      states{index} = {S L statesLS{i}{3:4} iJ 1 i};
  end
end

num_states = length(states);

for iq = 1:(2*k+1); Ukq{iq} = zeros(num_states); end;
%for i = 1:num_states
for i = 1:num_states
  L = states{i}{2}; S = states{i}{1}; J = states{i}{5}; Jz = states{i}{6}; irm = states{i}{7};
  %for j = 1:num_states
  for j = i:num_states   % Calculate only upper triangle. Time per q-value: 533.56s
    Lp = states{j}{2}; Sp = states{j}{1}; Jp = states{j}{5}; Jzp = states{j}{6}; irmp = states{j}{7};
    for iq = indq
      q = iq-1-k;
%      rm = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]) * redmat(irm,irmp);
      % NB. wigner symbol (j1j2jm|j1j2m1m2) = (-1)^(j1-j2-m) * sqrt(2*j+1) * [j1 j2 j; m1 m2 m]
      %Ukq{iq}(i,j) = (-1)^(J+Jz+k+q) * rm * threej([J Jp -Jz; k -q Jzp]) / sqrt(2*k+1);
      if S==Sp
        % Using Eqn 4 of Chan and Lam, 1970
%        Ukq{iq}(i,j) = (-1)^(k+S+Lp+(2*J)-Jz) * sqrt((2*J+1)*(2*Jp+1)) * threej([J k Jp; -Jz q Jzp]) ...
%	               * sixj([L J S; Jp Lp k]) * redmat(irm,irmp);
        rm = (-1)^(k+S+Lp+J) * sqrt((2*J+1)*(2*Jp+1)) * sixj([L J S; Jp Lp k]) * redmat(irm,irmp);
        % NB. The 3j symbol as I've program takes as arguments: [j1 j2 j3; J1 J2 -J3]!!
        % So in notations, the 3j below should be: 3j([J k Jp; -Jz q Jzp])
%        if mod(n,2)==1
	  %rm2 = real( (-1)^(k+S+Lp+2*J-Jz) * sqrt((2*J+1)*(2*Jp+1)) * sixj([L J S; Jp Lp k]) * redmat(irm,irmp) ...
	  %            * threej([J k Jp; -Jz q -Jzp]) );
	  rm2 = real( (-1)^(J-Jz) * rm * threej([J k Jp; -Jz q -Jzp]) );
          Ukq{iq}(i,j) = rm2;
	  %if J~=Jp
          %  %Ukq{iq}(i,j) = min([sign(Jz) sign(Jzp)])*rm2;
          %  Ukq{iq}(i,j) = sign(Jz)*rm2;
          %else
          %  Ukq{iq}(i,j) = sign(Jz)*sign(Jzp)*rm2;
          %end
	  %%if Jz>=0 & (sign(Jz)==sign(Jzp))
          %%  Ukq{iq}(i,j) = rm2;
          %%else
          %%  Ukq{iq}(i,j) = -rm2;
          %%end
%        else
          %Ukq{iq}(i,j) = real( (-1)^(2*J-Jz-Jzp) * rm * threej([J k Jp; -Jz q -Jzp]) );
%          Ukq{iq}(i,j) = real( (-1)^(2*J-Jz-Jzp) * rm * threej([J k Jp; -Jz q -Jzp]) );
%        end
        %if mod(k,2)==1 & (Jz<0 & q==-1) | (Jz>0 & q==1) | (Jz<0 & q==0)
        %  Ukq{iq}(i,j) = -Ukq{iq}(i,j);
        %end
      end
    end
  end
end

for iq = indq
  Ukq{iq} = Ukq{iq} + Ukq{iq}' - eye(size(Ukq{iq})).*Ukq{iq};
end
