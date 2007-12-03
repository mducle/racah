function Vkkk = racah_Vkkk(n,l,K,states)

% The matrix elements for the three particle operators V(k,k',k'') is:
%
%          (k) (k'') (k') (k') (0)          ---           { k  k' k'' }
% (psi||({u   u     }    u    )   ||psi') = >   ___ ____  { _  _      }
%          h   i          j                 --- psi,psi'  { L' L  L   }
%                                         
%                          (k)  ___   ___   (k'')  ____   ____   (k')
%               x   (psi||u   ||psi) (psi||u     ||psi') (psi'||u    ||psi')
%
%                     { J1 J2 J3 }
% Now, the 6j symbols { L1 L2 L3 } must have for the triads (J1 J2 J3), (J1 L2 L3), (L1 J2 L3), (L1 L2 J3)
%   1. All satisfy the triangular inequality: |x-y| <= z <= x+y  
%   2. Elements of each triad sum to an integer
% So for non-zero elements we need:                 _        _
%   1. The 6j symbol to be nonzero, so the triads (kLL) and (L'k'L) to be nonzero. Thus we want:
%      ( |k-Lbar| < L < k+Lbar ) and ( |L'bar-k'| < L < L'bar+k )

if ~exist('states')
  states = racah_states(n,l);
end

num_states = length(states); Vkkk = zeros(num_states);
k   = K(1);             Uk   = racah_Umat(n,l,k,states);                         Vk = sqrt(2*k+1).*Uk;
kp  = K(2); if kp~=k;   Ukp  = racah_Umat(n,l,kp,states);  else Ukp = Uk; end;   Vkp = sqrt(2*kp+1).*Ukp;
kpp = K(3); if kpp~=kp; Ukpp = racah_Umat(n,l,kpp,states); else Ukpp = Ukp; end; Vkpp = sqrt(2*kpp+1).*Ukpp;

for i = 1:num_states
  L(i) = racah_lconv(states{i}{2});
end 

for i = 1:num_states
  for j = 1:num_states
    for b = 1:num_states   % Loop over psi_bar                             _
      %if (L(i) < (k+L(b))) | (L(i) > abs(k-L(b)))   % Ensures the triad (kLL) is non-zero
      for p = 1:num_states % Loop over psi'_bar                             _
      %if (L(i) < (kp+L(p))) | (L(i) > abs(kp-L(p))) % Ensures the triad (k'L'L) is non-zero
        Vkkk(i,j) = Vkkk(i,j) + ( sixj([k kp kpp; L(p) L(b) L(i)]) * Uk(i,b)*Ukpp(b,p)*Ukp(p,j) );
      %end % if (k'L'L)
      end % for p
      %end % if (kLL)
    end  % for b
    Uk2   = (1/(2*L(i)+1)) * ( (-1).^(L(i)-L) .* Uk(i,:)*Uk(:,i) ); 
    Ukp2  = (1/(2*L(i)+1)) * ( (-1).^(L(i)-L) .* Ukp(i,:)*Ukp(:,i) ); 
    Ukpp2 = (1/(2*L(i)+1)) * ( (-1).^(L(i)-L) .* Ukpp(i,:)*Ukpp(:,i) ); 
    if (i==j); delta = (2*n)/(2*l+1); else delta = 0; end
    Vkkk(i,j) = Vkkk(i,j) - ( sixj([k kp kpp; l l l]) * ( Uk2 + Ukp2 + Ukpp2 - delta ) );
  end  % for j
end  % for i

Vkkk = Vkkk .* ( sqrt(2*k+1) * sqrt(2*kp+1) );
