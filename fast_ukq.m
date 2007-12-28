function [Ukq,states] = fast_ukq(n,l,k,indq)
% Calculates the matrix Ukq, after Elliot, Judd, and Runciman using a fast algorithm
% 
% Syntax: [Ukq,states] = fast_ukq(n,l,k,q)
%
% Inputs:  n - scalar - number of equivalent "l"-electrons
%          l - scalar - orbital angular momentum number of each equivalent electron
%          k - scalar - rank of Ukq tensor (only 2, 4, 6 are valid).
%          q - scalar - (optional) the order of the component of Ukq tensor. If not
%                       specified, a cell array will be returned with all the 2k+1 
%                       components -k,-k+1,...,k.
%
% Outputs: Ukq    - matrix  - the Ukq sparse matrix.
%          states - cell    - cell with the [S L U v J mJ] values of each state

% By Duc Le - 2007 - duc.le@ucl.ac.uk

% The function calculates first the reduce matrix elements <vULS||U^k||v'U'L'S>
% using the racah_Umat function. It then only calculates the subsequent J- and mJ-
% dependent matrix elements for the values of the reduced matrix which is non-zero.
% It does this by calculating a J-dependent only matrix <vULSJ||U^k||v'U'L'SJ'> 
% for each non-zero reduced matrix element, and then an mJ-dependent matrix 
% <vULSJmJ|Ukq|v'U'L'SJ'mJ> for each J-J' value.
% 
% The function then slots each of these smaller "J-blocks" into the full Ukq matrix
% leaving the rest zeros.
%
% Time to calculate a single q-value matrix is approx 0.3s for f2 and approx 3min
% for f5 on a Pentium M 1.6GHz. The original racah_Ukq function which calculated
% every matrix element took approx 6s for f2 and approx 45min for f5. 

if ~exist('indq')
  indq = 1:(2*k+1);
else
  indq = indq + k+1;
end

statesLS = racah_states(n,l);
redmat = racah_Umat(n,l,k,statesLS);

% Determines the L S J values for each matrix elements and the index of each J-J' block
num_states = 0; LStest = [0 0];
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2}); Lnum(i) = L;
  S = statesLS{i}{1}; Snum(i) = S;
  Jmin = abs(L-S); Jmax = L+S;
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
end

% Only calculate if the user wants the list of states.
if (nargin>1)
  index = 0;
  for i = 1:length(statesLS)
    for iJ = Jmin:Jmax
      for imJ = -Jmin:Jmax
        index = index + 1;
        states{index} = {S L statesLS{i}{3:4} iJ imJ i};
      end
    end
  end
end

% Determines the J-values for which we need to calculate operator equivalent matrices
minJ = min(Jcount(:,2)); maxJ = max(Jcount(:,3)); valJ = minJ:maxJ;

% Initialises the Ukq matrix to zeros using a sparse matrix for smaller storage - important for n>4
if size(indq~=1)
  for iq = indq; Ukq{iq} = sparse(num_states,num_states); end;
else
  Ukq{iq} = sparse(num_states,num_states);
end

for iq = indq
  q = iq-1-k;

  % Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
  %   each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
  for iJ = 1:length(valJ)
    for iJp = 1:length(valJ)
      mJmat{iJ,iJp} = 0;
    end
  end

  % Only calculates for the non-zero values of the reduced matrix elements <vULS||U^k||v'U'L'S>
  [rm_i,rm_j] = find(redmat);
  for irm = 1:length(rm_i)
    rmi = rm_i(irm); rmj = rm_j(irm);
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

      Ukq{iq}(indexJstart(rmi):indexJstop(rmi),indexJstart(rmj):indexJstop(rmj)) = rmJ .* redmat(rmi,rmj);

    end  % if S==Sp
  end    % for irm
end      % for iq

if length(indq)==1
  Ukq = Ukq{indq};
end
