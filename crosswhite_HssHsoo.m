function Hss = crosswhite_HssHsoo(n,l,M,P)

if l~=3; error('Sorry, only l=3 implemented'); end

if length(M)==6
  MP = M;
elseif length(M)==3 || length(P)==3
  MP = [M P];
else
  error('M must be a length 3 or length 6 vector; P must be a length 3 vector');
end

switch n
  case 2; crosswhite_f2MPmat; redmat = {M0mat M2mat M4mat P2mat P4mat P6mat};
  case 3; crosswhite_f3MPmat; redmat = {M0mat M2mat M4mat P2mat P4mat P6mat};
  case 4; crosswhite_f4MPmat; redmat = {M0mat M2mat M4mat P2mat P4mat P6mat};
  case 5; crosswhite_f5MPmat; redmat = {M0mat M2mat M4mat P2mat P4mat P6mat};
otherwise
  error('Sorry, only n=2-5, l=3 implemented so far');
end

statesLS = racah_states(n,l); index = 0;
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2});
  S = statesLS{i}{1};
  Jmin = abs(L-S); Jmax = L+S;
  for iJ = Jmin:Jmax
    index = index + 1;
    states{index} = {S statesLS{i}{2:4} iJ i};
  end
  LSdiff(i) = L-S; 
end
for i = 1:length(statesLS); lTriFact(:,i) = LSdiff(i)-LSdiff; end; lTriFact = tril((-1).^lTriFact,-1);
for i = 1:6; redmat{i} = triu(redmat{i}) + triu(redmat{i})'.*lTriFact; end

% Equation for the matrix elements is:
%                                      S'+L+J { S' L' J }         (tt)
% <vUSLJ|H|v'U'S'L'J'> = d(J,J') x (-1)       { L  S  t } <vUSL||T    ||v'U'S'L'>
%
% So, to satisfy the six-j symbols we have to have the triads (S'St) and (LL't) non-zero
% i.e. |S-S'|<= t <= S+S' and |L-L'|<= t <= L+L' with t being either 1 or 2.
%
% Ref: Judd, Crosswhite and Crosswhite, Phys. Rev. 169(1968)130
num_states = length(states); Hss = sparse(num_states,num_states);
for i = 1:num_states
  S = states{i}{1}; Ls = states{i}{2}; L = racah_lconv(Ls); J = states{i}{5}; irm = states{i}{6};
  for j = 1:num_states
    Sp = states{j}{1}; Lsp = states{j}{2}; Lp = racah_lconv(Lsp); Jp = states{j}{5}; jrm = states{j}{6};
    if J==Jp %&& (S-Sp)<=2 && (L-Lp)<=2
      %for z = 1:3
      %  mk{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 2]) * redmat{z}(irm,jrm);
      %end
      if (S-Sp)<=1 && (L-Lp)<=1
        %for z = 4:6
        for z = 1:6
          Hss(i,j) = Hss(i,j) + MP(z) * ( (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 1]) * redmat{z}(irm,jrm) );
        end
        %for z = 1:4
        %  pk{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 1]) * redmat{z+6}(irm,jrm);
        %end
      end 
    end % if J==Jp
  end
end
