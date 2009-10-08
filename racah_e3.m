function e3 = racah_e3(n,l,states)
% Calculates the Racah electrostatic energy operator e3, after Racah IV.

if ~exist('states')
  states = racah_states(n,l);
end

if n>7; n=14-n; end

for j = 1:length(states)
  for i = 1:j %length(states)

    Si = states{i}{1}; Li = states{i}{2}; vi = states{i}{3}; Ui = states{i}{4};
    Sj = states{j}{1}; Lj = states{j}{2}; vj = states{j}{3}; Uj = states{j}{4};
    Wi = racah_vtow(vi,Si,l); Wj = racah_vtow(vj,Sj,l);
    
    if Si==(n/2) & i==j               % States of maximum spin, diagonal elements
      L = racah_lconv(Li);
      e3(i,i) = -3 * (L*(L+1)/2 - 12*racah_g(Ui));   % Eqn 81 of Racah IV
    %elseif Li==Lj & Si==Sj
    elseif racah_lconv(Li)==racah_lconv(Lj) & Si==Sj
      if vi == vj
        v = vi;
        if (n==6 & v==6) | (n==7 & v==7)             % Eqn 83 of Racah IV      
          e3(i,j) = 0;
        elseif n==(v+2)                              % Eqn 82 and 84a of Racah IV
          e3(i,j) = ((1-v)/(7-v)) * ( racah_yfn(vi,vi,Si,Ui,vi,Uj) * racah_phi(Ui,Uj,Li,Lj) );
        elseif n==(v+4)                              % Eqn 82 and 84b of Racah IV
          e3(i,j) = (-4/(7-v))    * ( racah_yfn(vi,vi,Si,Ui,vi,Uj) * racah_phi(Ui,Uj,Li,Lj) );
        else                                         % Eqn 87 of Racah IV
          e3(i,j) = racah_yfn(n,v,Si,Ui,v,Uj) * racah_phi(Ui,Uj,Li,Lj);
	end
      elseif n==5 & vi==1 & vj==3 & Si==1/2          % Eqn 86a of Racah IV
        e3(i,j) = sqrt(2/5) * ( racah_yfn(3,1,1/2,Ui,3,Uj) * racah_phi(Ui,Uj,Li,Lj) );
      elseif n==6 & vi==0 & vj==4 & Si==0            % Eqn 86b of Racah IV
        e3(i,j) = sqrt(9/5) * ( racah_yfn(4,0,0,Ui,4,Uj)   * racah_phi(Ui,Uj,Li,Lj) );
      elseif n==6 & vi==2 & vj==4                    % Eqn 86c of Racah IV
        e3(i,j) = sqrt(1/6) * ( racah_yfn(4,2,Si,Ui,4,Uj)  * racah_phi(Ui,Uj,Li,Lj) );
      elseif n==7 & vi==1 & vj==5 & Si==1/2          % Eqn 86d of Racah IV
        e3(i,j) = sqrt(3/2) * ( racah_yfn(5,1,1/2,Ui,5,Uj) * racah_phi(Ui,Uj,Li,Lj) );
      else                                           % Eqn 87 of Racah IV
        e3(i,j) = racah_yfn(n,vi,Si,Ui,vj,Uj) * racah_phi(Ui,Uj,Li,Lj);
      end
      if i==j && e3(i,i)~=0
        L = racah_lconv(Li);
        e3(i,i) = e3(i,i) - (L*(L+1)/2 - 12*racah_g(Ui));
      end
    else
      e3(i,j) = 0;
    end
  end
end

% Adds in the lower triangle to make e3 hermitian
e3diag = eye(size(e3)).*e3;
e3 = e3 + e3' - e3diag;
