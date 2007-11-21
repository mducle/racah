function [Ukq,states] = racah_Ukq(n,l,k,indq)
% Calculates the matrix Ukq, after Elliot, Judd, and Runciman.

icfact = [-sqrt(15/7)/2       sqrt(11/14)        -sqrt(429/7)/10  ];

if ~exist('indq')
  indq = 1:(2*k+1);
else
  indq = indq + k+1;
end

statesLS = racah_states(n,l);
redmat = racah_Umat(n,l,k,statesLS) ./ icfact(k/2);

%indLS = 0; index = 0; 
num_states = 0; LStest = [0 0];
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2}); Lnum(i) = L;
  S = statesLS{i}{1}; Snum(i) = S;
  Jmin = abs(L-S); Jmax = L+S; %Jind = Jmin:Jmax;
  Jcount(i,:) = [Jmax-Jmin Jmin Jmax];
  indexJstart(i) = num_states + 1;
  num_states = num_states + sum(2.*[Jmin:Jmax]+1);    % Total number of |vUSLJmJ> states
  indexJstop(i) = num_states;
  if [L S]~=LStest
    uniqueLS(i) = 0;                                  % So we know which values of LSJ
  else                                                %   to calculate the reduced
    uniqueLS(i) = 1;                                  %   matrix elements <J||U^k||J'>
  end                                                 %   for (see below).
  LStest = [L S];
    
  %indexJstart(i) = index; index = index+Jcount(i,1)-1; indexJstop(i) = index;
  %for iJ = 1:Jcount(i,1)
  %  indexmJstart{i}(iJ) = numstates;                  % numstates is running total
  %  num_states = num_states + (2*(Jmin+iJ-1)+1);      %   number of |vUSLJmJ> states
  %  indexmJstop{i}(iJ) = numstates;
  %end

  %for iJ = Jmin:Jmax
  %  for imJ = -Jmin:Jmax
  %    index = index + 1;
  %    states{index} = {S L statesLS{i}{3:4} iJ imJ i};
  %  end
  %end
end

% Determines the J-values for which we need to calculate operator equivalent matrices
minJ = min(Jcount(:,2)); maxJ = max(Jcount(:,3)); valJ = minJ:maxJ;

% Initialises the Ukq matrix to zeros
if size(indq~=1)
  %for iq = indq; Ukq{iq} = zeros(num_states); end;
  for iq = indq; Ukq{iq} = sparse(num_states,num_states); end;
else
  %Ukq{iq} = zeros(num_states);
  Ukq{iq} = sparse(num_states,num_states);
end

for iq = indq
  q = iq-1-k;
  % Uses operator equivalents to calculate matrix elements of mJ values within a J-manifold.
  %for iJ = 1:length(valJ)
  %  J = valJ(iJ);
  %  for imJ_i = 1:(2*J+1)
  %    Jz = imJ_i-J-1;
  %    for imJ_j = 1:(2*J+1)
  %      Jzp = imJ_j-J-1;
  %      mJmat{iJ}(imJ_i,imJ_j) = (-1)^(J+Jz+k+q) * (1/sqrt(2*k+1)) * wigner([J Jp -Jz Jzp J Jp k -q]);
  %    end
  %  end
  %end

  % Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
  %   each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
  for iJ = 1:length(valJ)
    for iJp = 1:length(valJ)
      mJmat{iJ,iJp} = 0; %zeros(2*valJ(iJ)+1,2*valJ*iJp+1);
    end
  end

  % Caculate the J-dependent matrix elements <J||U^k||J'> = <vULSJ||U^k||v'U'L'SJ'>/<vULS||U^k||v'U'L'S>
  %index = 0;
  %for iLS = 1:size(uniqueLS,1)
  %  L = unique(iLS,1); S = unique(iLS,2);
  %  Jmin = uniqueLS(iLS,3); Jmax = uniqueLS(iLS,4); Jval = [Jmin:Jmax];
  %  for iJ_i = 1:length(Jval)
  %    J = Jval(iJ_i);
  %    for iJ_j = 1:length(Jval)
  %      Jp = Jval(iJ_i);
  %      Jmat{index}(iJ_i,iJ_j) = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]);
  %    end
  %  end
  %  if unique(iLS,5)==0
  %    index = index + 1;
  %  else
  %    for i = 1:unique(iLS,5)                         % To avoid repeatedly calculating the same
  %      Jmat{index+i} = Jmat{index};                  %   reduced matrix for save LSJ's in cases
  %    end                                             %   where different states have same LSJ but
  %    index = index + 1 + unique(iLS,5);              %   different vU's.
  %  end
  %end

  % Only calculates for the non-zero values of the reduced matrix elements <vULS||U^k||v'U'L'S>
  [rm_i,rm_j] = find(redmat);
  for irm = 1:length(rm_i)
    rmi = rm_i(irm); rmj = rm_j(irm);
    %L = Lnum(rm_i(irm)); S = Snum(rm_i(irm)); Lp = Lnum(rm_j(irm)); Sp = Snum(rm_j(irm));
    L = Lnum(rmi); S = Snum(rmi); Lp = Lnum(rmj); Sp = Snum(rmj);
    if S==Sp

      % Caculate the J-dependent matrix elements <J||U^k||J'> = <vULSJ||U^k||v'U'L'SJ'>/<vULS||U^k||v'U'L'S>
      Jmin  = abs(L-S);  Jmax  = L+S;  Jval  = [Jmin:Jmax];   count = 0;
      Jminp = abs(Lp-S); Jmaxp = Lp+S; Jvalp = [Jminp:Jmaxp]; countp = 0;
      rmJ = [];
      for iJ_i = 1:length(Jval)
        J = Jval(iJ_i); countp = 0;
        for iJ_j = 1:length(Jvalp)
	  Jp = Jvalp(iJ_j);
	  rm = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]);

          % Calculates the q- and Jz- dependent matrix elements <mJ|Ukq|mJ'> = <psi|Ukq|psi'>/<J||U^k||J'>
	  valJ_i = J-minJ+1; valJ_j = Jp-minJ+1;
	  if mJmat{valJ_i,valJ_j} == 0
	    for imJ_i = 1:(2*J+1)
	      Jz = imJ_i-1-J;
	      for imJ_j = 1:(2*Jp+1)
	        Jzp = imJ_j-1-Jp;
                mJmat{valJ_i,valJ_j}(imJ_i,imJ_j) = (-1)^(J+Jz+k+q) ...
		                                    * (1/sqrt(2*k+1)) * wigner([J Jp -Jz Jzp J Jp k -q]);
              end
	    end
	  end
	  rmJ((count+1):(count+2*J+1),(countp+1):(countp+2*Jp+1)) = mJmat{valJ_i,valJ_j} .* rm;
          countp = countp + (2*Jp+1);
        end
        count = count + (2*J+1);
      end



      %---
      %if unique(iLS,5)==0
      %  index = index + 1;
      %else
      %  for i = 1:unique(iLS,5)                       % To avoid repeatedly calculating the same
      %    Jmat{index+i} = Jmat{index};                %   reduced matrix for save LSJ's in cases
      %  end                                           %   where different states have same LSJ but
      %  index = index + 1 + unique(iLS,5);            %   different vU's.
      %end

      %Ukq{iq}(indexJstart(rm_i(irm)):indexJstop(rm_i(irm)),indexJstart(rm_j(irm)):indexJstop(rm_j(irm))) = ...
      %         rmJ .* redmat(rm_i(irm),rm_j(irm));

      Ukq{iq}(indexJstart(rmi):indexJstop(rmi),indexJstart(rmj):indexJstop(rmj)) = rmJ .* redmat(rmi,rmj);

    end  % if S==Sp
  end    % for irm
  
  %Ukq{iq} = Ukq{iq} ./ icfact(k/2);

end      % for iq
    

%----------------

%num_states = length(states);
%
%for i = 1:num_states
%  L = states{i}{2}; S = states{i}{1}; J = states{i}{5}; Jz = states{i}{6}; irm = states{i}{7};
%  for j = 1:num_states 
%    Lp = states{j}{2}; Sp = states{j}{1}; Jp = states{j}{5}; Jzp = states{j}{6}; irmp = states{j}{7};
%    rm = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]) * redmat(irm,irmp);
%    for iq = indq
%      q = iq-1-k;
%      % NB. wigner symbol (j1j2jm|j1j2m1m2) = (-1)^(j1-j2-m) * sqrt(2*j+1) * [j1 j2 j; m1 m2 m]
%      if S==Sp
%        Ukq{iq}(i,j) = (-1)^(J+Jz+k+q) * (1/sqrt(2*k+1)) * rm * wigner([J Jp -Jz Jzp J Jp k -q]) / icfact(k/2);
%      end
%    end
%  end
%end

if length(indq)==1
  Ukq = Ukq{indq};
end
