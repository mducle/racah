function [gs,es] = hundsrule(n,l,flag)
% Calculates the Hund's rule ground state and excited state Spin-Orbit terms.

% Finds the ground state
% Rule 1: Maximise S (or multiplicity 2S+1)
if n<(2*l+1)
  S = n/2;            % Less than half full - all spins add.
else
  S = (n-(2*l+1))/2;  % More than half full - some spins cancel.
end

% Rule 2: Maximise L
ml = l:-1:-l;
L = sum(ml(1:n)); 

% Rule 3: Minimise J if shell less than half full; Maximise J otherwise
if n<(2*l+1)
  J = abs(L-S);       % Less than half full.
else
  J = L+S;            % More than half full.
end

gs = [L S J];

% Outputs Terms Symbol instead
if nargin==3
  switch L
    case 0; gs = sprintf('^{%1.0f}S_{%1.1f}',[2*S+1 J]);
    case 1; gs = sprintf('^{%1.0f}P_{%1.1f}',[2*S+1 J]);
    case 2; gs = sprintf('^{%1.0f}D_{%1.1f}',[2*S+1 J]);
    case 3; gs = sprintf('^{%1.0f}F_{%1.1f}',[2*S+1 J]);
    case 4; gs = sprintf('^{%1.0f}G_{%1.1f}',[2*S+1 J]);
    case 5; gs = sprintf('^{%1.0f}H_{%1.1f}',[2*S+1 J]);
    case 6; gs = sprintf('^{%1.0f}I_{%1.1f}',[2*S+1 J]);
    case 7; gs = sprintf('^{%1.0f}J_{%1.1f}',[2*S+1 J]);
    case 8; gs = sprintf('^{%1.0f}K_{%1.1f}',[2*S+1 J]);
    case 9; gs = sprintf('^{%1.0f}L_{%1.1f}',[2*S+1 J]);
  end
end

if nargout == 2

% List the excited state in LS-coupling order.

Smax = S;
if mod(n,2)==1
  Smin = 1/2;          % Odd number of electrons.
else
  Smin = 0;            % Even number of electrons.
end

index = 1;
S = Smax:-1:Smin;
for iS = 1:length(S)

  n1 = n-(Smax-S(iS));
  Lmax = sum(ml(1:n1)) + sum(ml(1:(n-n1)));
  if mod(n,2)==1; Lmin = 1; else; Lmin = 0; end
  L = Lmax:-1:Lmin;
  for iL = 1:length(L)
  
    if n<(2*l+1)
      J = abs(L(iL)-S(iS)):(L(iL)+S(iS));
    else 
      J = (L(iL)+S(iS)):-1:abs(L(iL)-S(iS))
    end
    for iJ = 1:length(J)
      es{index} = [L(iL) S(iS) J(iJ)];                % Lists be increasing J values within Term.
      index = index+1;
    end
  
  end
end

if nargin==3
  for in = 1:(index-1)
    switch es{in}(1)
      case 0; es2{in} = sprintf('^{%1.0f}S_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 1; es2{in} = sprintf('^{%1.0f}P_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 2; es2{in} = sprintf('^{%1.0f}D_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 3; es2{in} = sprintf('^{%1.0f}F_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 4; es2{in} = sprintf('^{%1.0f}G_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 5; es2{in} = sprintf('^{%1.0f}H_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 6; es2{in} = sprintf('^{%1.0f}I_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 7; es2{in} = sprintf('^{%1.0f}J_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 8; es2{in} = sprintf('^{%1.0f}K_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
      case 9; es2{in} = sprintf('^{%1.0f}L_{%1.1f}', [2*es{in}(2)+1 es{in}(3)]);
    end
  end
  es = es2;
end

end
