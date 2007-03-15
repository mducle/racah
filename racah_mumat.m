function mu = racah_mumat(n,l)
% Calculates the magnetic moment operator matrices mu_x, mu_y, mu_z in the |vSLJM> basis

% Constants
g_s = 2.0;     % electronic gyromagnetic ratio

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

for iq = 1:3; Ukq{iq} = zeros(num_states); end;
for i = 1:num_states
  L = states{i}{2}; S = states{i}{1}; J = states{i}{5}; M = states{i}{6};
  if J>1/2
    g = 1 + (g_s-1) * (J*(J+1) - L*(L+1) + S*(S+1)) / (2*J*(J+1));
    f = sqrt( (S+L+J+1)*(S+L-J+1)*(S+J-L)*(L+J-S) / (J^2*(2*J+1)*(2*J-1)) );
    for j = 1:num_states
      Lp = states{j}{2}; Sp = states{j}{1}; Jp = states{j}{5}; Mp = states{j}{6}; 
      
      if L==Lp & S==Sp
        if J==Jp
          if M==Mp
            mu{3}(i,j) = M*g;                         % mu_z
          elseif M==(Mp+1)
            mu{1}(i,j) = sqrt( (J+M+1)*(J-M) ) * g/2; % mu_x 
            mu{2}(i,j) = sqrt(-(J+M+1)*(J-M) ) * g/2; % mu_y
          elseif M==(Mp-1)
            mu{1}(i,j) = sqrt( (J-M+1)*(J+M) ) * g/2;
            mu{2}(i,j) =-sqrt(-(J-M+1)*(J+M) ) * g/2;
          end
        elseif J==(Jp-1)
          if M==Mp
            mu{3}(i,j) = (g_s-1) * sqrt(J^2-M^2) * f/2;
          elseif M==(Mp+1)
            mu{1}(i,j) = (g_s-1) * sqrt( (J-M-1)*(J-M) ) * f/4;
            mu{2}(i,j) = (g_s-1) * sqrt( (J-M-1)*(J-M) ) * f/4;
          elseif M==(Mp-1)
            mu{1}(i,j) =-(g_s-1) * sqrt( (J+M-1)*(J+M) ) * f/4;
            mu{2}(i,j) = (g_s-1) * sqrt(-(J+M-1)*(J+M) ) * f/4;
          end
        elseif J==(Jp+1)
          if M==Mp
            mu{3}(i,j) = (g_s-1) * sqrt( (J+M+1)*(J-M+1) ) * f/2;
          elseif M==(Mp+1)
            mu{1}(i,j) =-(g_s-1) * sqrt( (J+M+1)*(J+M+2) ) * f/4;
            mu{2}(i,j) =-(g_s-1) * sqrt(-(J+M+1)*(J-M+2) ) * f/4;
          elseif M==(Mp-1)
            mu{1}(i,j) = (g_s-1) * sqrt( (J-M+1)*(J+M+2) ) * f/4;
            mu{2}(i,j) =-(g_s-1) * sqrt(-(J-M+1)*(J+M+2) ) * f/4;
          end
        end
      end

    end
  end
end

% Time to execute: n=5, t=535.2s
