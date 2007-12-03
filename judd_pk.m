function pk = judd_pk(n,l,Pk,cellflag)
% Calculates the spin-spin Hamiltonian for the l^n configuration by a chain-calculation

if ~isscalar(n) | ~isscalar(l) | ~isnumeric(n) | ~isnumeric(l) 
  error('n, l must be numerical scalars');
elseif ~isvector(M) | length(M)~=3
  error('M must be a three-component vector');
end

%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
px = [  1            0              0           1  sqrt(15/2)           0           0;      % 3P
        0     sqrt(14)              0           0    sqrt(10)    sqrt(11)           0;      % 3F
        0            0       sqrt(55)           0           0    sqrt(10)  sqrt(13/2);      % 3H
        0            0              0           0           0           0           0;      % 1S
        0            0              0           0           0           0           0;      % 1D
        0            0              0           0           0           0           0;      % 1G
        0            0              0           0           0           0           0];     % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
p0 = [ -1            0              0          -2           1           0           0;      % 3P
        0           -1              0           0          -1           1           0;      % 3F
        0            0             -1           0           0          -1           1;      % 3H
        0            0              0           0           0           0           0;      % 1S
        0            0              0           0           0           0           0;      % 1D
        0            0              0           0           0           0           0;      % 1G
        0            0              0           0           0           0           0].*mx; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
p2 = [-45            0              0        -105          32           0           0;      % 3P
        0           10              0           0        -9/2         -20           0;      % 3F
        0            0             25           0           0        55/2           0;      % 3H
        0            0              0           0           0           0           0;      % 1S
        0            0              0           0           0           0           0;      % 1D
        0            0              0           0           0           0           0;      % 1G
        0            0              0           0           0           0           0].*mx; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
p4 = [-33            0              0        -231         -33           0           0;      % 3P
        0           33              0           0          66          32           0;      % 3F
        0            0             51           0           0         -23         -21;      % 3H
        0            0              0           0           0           0           0;      % 1S
        0            0              0           0           0           0           0;      % 1D
        0            0              0           0           0           0           0;      % 1G
        0            0              0           0           0           0           0].*mx; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
p6 = [1287           0              0        -429        -286           0           0;      % 3P
        0          286              0           0      -429/2        -104           0;      % 3F
        0            0             13           0           0       -65/2          -6;      % 3H
        0            0              0           0           0           0           0;      % 1S
        0            0              0           0           0           0           0;      % 1D
        0            0              0           0           0           0           0;      % 1G
        0            0              0           0           0           0           0].*mx; % 1I

p0 = p0+triu(p0,1)'; p2 = p2+triu(p2,1)'; p4 = p4+triu(p4,1)'; p6 = p6+triu(p6,1)';

oldp = {p0 p2 p4 p6}; st_old = racah_states(2,l); s = 1/2; t = 1;

if n==2 
  Hss = p0.*Pk(1) + p2.*Pk(2) + p4.*Pk(3) + p6.*Pk(4);
  if nargin~=4
    return;
  else
    Hss.f2 = {p0 p2 p4 p6 Hss};
    return;
  end
elseif n<2 || n>14
  error('number of equivalent electron, n, must be between 2 and 14');
end

for n_calc = 3:n
    
  display(sprintf('Calculating matrix for %c^%1g',lower(racah_lconv(l)),n_calc)); tic

  st_now = racah_states(n_calc,l); num_states = length(st_now); clear cfp;
  for i = 1:num_states
    S(i) = st_now{i}{1}; Ls = st_now{i}{2}; v  = st_now{i}{3}; U  = st_now{i}{4};
    L(i) = racah_lconv(Ls); LSdiff(i) = L(i)-S(i);
    [coefs,parents,ind_parents]= racah_parents(n_calc,l,v,U,S(i),Ls);
    cfp{i} = zeros(1,length(st_old)); cfp{i}(ind_parents) = coefs;
  end
  clear lTriFact;
  for i = 1:num_states
    lTriFact(:,i) = LSdiff(i)-LSdiff;
  end
  lTriFact = tril((-1).^lTriFact,-1);
  for i = 1:length(st_old)
    Sb(i) = st_old{i}{1}; Lb(i) = racah_lconv(st_old{i}{2});
  end  

  for im = 1:4
    pk{im} = zeros(num_states);
    for i = 1:num_states        % Calculates only the upper triangle
      for j = i:num_states      % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
        if abs(S(j)-S(i))<=1 && abs(L(j)-L(i))<=1
        cfplogmat = cfp{i}'*cfp{j}; [inz,jnz] = find( cfplogmat.*oldm11{im} );
        for p = 1:length(inz);
          pi = inz(p); pj = jnz(p);
          mat_el = cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
                   * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1));
          pk{im}(i,j) = pk{im}(i,j) + ( mat_el * sixj([S(i) t S(j); Sb(pj) s Sb(pi)]) ...
                                               * sixj([L(i) t L(j); Lb(pj) l Lb(pi)]) * oldp{im}(pi,pj) );
        end
        end % if /\S=0,1, /\L=0,1
      end
    end
    pk{im} = pk{im}+pk{im}'.*lTriFact; pk{im} = ((n_calc)/(n_calc-2)) .* pk{im};
  end

  time = toc; if time<60; display(sprintf('Elapsed time was %0.2g s',time));
                     else display(sprintf('Elapsed time was %0.2g min',time/60)); end;

  %if exist('matname') && ~exist(matfilename,'file')
  %  display(sprintf('Saving %s',matfilename));
  %  save(matfilename,'matrix');
  %end

  if nargin==4
    H = pk{1}.*Pk(1) + pk{2}.*Pk(2) + pk{3}.*Pk(3) + pk{4}.*Pk(4);
    Hss.(sprintf('%c%1g',lower(racah_lconv(l)),n_calc)) = {pk{1} pk{2} pk{3} pk{4} H};
  end
  
  oldp = pk;
  st_old = st_now;

end

if nargin~=4
  Hss = pk{1}.*Pk(1) + pk{2}.*Pk(2) + pk{3}.*Pk(3) + pk{4}.*Pk(4);
end