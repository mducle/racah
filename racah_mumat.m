function mu = racah_mumat(n,l)
% Calculates the magnetic moment operator matrices mu_x, mu_y, mu_z in the |vSLJM> basis

% Physical constants. Taken from NIST Reference on Constants, Units, and 
% Uncertainty, http://physics.nist.gov/cuu/Constants/
g_s = 2.0023193043622;     % electronic g-factor

statesLS = racah_states(n,l);

index = 0;
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2});
  S = statesLS{i}{1};
  Jmin = abs(L-S); Jmax = L+S;
  for iJ = Jmin:Jmax
    for mJ = -iJ:iJ
      index = index + 1;
      states{index} = {S L statesLS{i}{3:4} iJ mJ};
    end
  end
end

num_states = length(states);

for iq = 1:3; mu{iq} = zeros(num_states); end;
for i = 1:num_states
  L = states{i}{2}; S = states{i}{1}; J = states{i}{5}; M = states{i}{6};
  g = 1 + (g_s-1) * (J*(J+1) - L*(L+1) + S*(S+1)) / (2*J*(J+1));
  denom = (J^2*(2*J+1)*(2*J-1));
  if denom~=0
    f  = sqrt( (S+L+J+1)*(S+L-J+1)*(S+J-L)*(L+J-S) / denom );
%   f  = sqrt( (S+L+J+1)*(S+L-J+1)*(S+J-L)*(L+J-S) / (J^2*(2*J+1)*(2*J-1)) );
  else
    f = 0;
  end
  denom = ((J+1)^2*(2*(J+1)+1)*(2*(J+1)-1));
  if denom~=0
    fp = sqrt( (S+L+(J+1)+1)*(S+L-(J+1)+1)*(S+(J+1)-L)*(L+(J+1)-S) / denom );
%   fp = sqrt( (S+L+(J+1)+1)*(S+L-(J+1)+1)*(S+(J+1)-L)*(L+(J+1)-S) / ((J+1)^2*(2*(J+1)+1)*(2*(J+1)-1)) );
  else
    fp = 0;
  end
  for j = 1:num_states
    Lp = states{j}{2}; Sp = states{j}{1}; Jp = states{j}{5}; Mp = states{j}{6}; 
      
    if states{i}{3}==states{j}{3} & states{i}{4}==states{j}{4} & L==Lp & S==Sp
      if J==Jp
        if M==Mp
          mu{3}(i,j) = M*g;                         % mu_z
        elseif Mp==(M+1)
          mu{1}(i,j) = sqrt( (J+M+1)*(J-M) ) * g/2; % mu_x 
          mu{2}(i,j) = sqrt(-(J+M+1)*(J-M) ) * g/2; % mu_y
        elseif Mp==(M-1)
          mu{1}(i,j) = sqrt( (J-M+1)*(J+M) ) * g/2;
          mu{2}(i,j) =-sqrt(-(J-M+1)*(J+M) ) * g/2;
        end
      elseif Jp==(J-1)
        if M==Mp
          mu{3}(i,j) = (g_s-1) * sqrt(J^2-M^2) * f/2;
        elseif Mp==(M+1)
          mu{1}(i,j) = (g_s-1) * sqrt( (J-M-1)*(J-M) ) * f/4;
          mu{2}(i,j) = (g_s-1) * sqrt(-(J-M-1)*(J-M) ) * f/4;
        elseif Mp==(M-1)
          mu{1}(i,j) =-(g_s-1) * sqrt( (J+M-1)*(J+M) ) * f/4;
          mu{2}(i,j) = (g_s-1) * sqrt(-(J+M-1)*(J+M) ) * f/4;
        end
      elseif Jp==(J+1)
        if M==Mp
          mu{3}(i,j) = (g_s-1) * sqrt( (J+M+1)*(J-M+1) ) * fp/2;
        elseif Mp==(M+1)
          mu{1}(i,j) =-(g_s-1) * sqrt( (J+M+1)*(J+M+2) ) * fp/4;
          mu{2}(i,j) =-(g_s-1) * sqrt(-(J+M+1)*(J+M+2) ) * fp/4;
        elseif Mp==(M-1)
          mu{1}(i,j) = (g_s-1) * sqrt( (J-M+1)*(J-M+2) ) * fp/4;
          mu{2}(i,j) =-(g_s-1) * sqrt(-(J-M+1)*(J-M+2) ) * fp/4;
        end
      end
    end

  end
end

% Time to execute: n=5, t=535.2s
