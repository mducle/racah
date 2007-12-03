function T22 = judd_T22(n,l)
% Calculates the z_i operators (i=1-4) that forms the T^(22) tensor representing the spin-spin Hamiltonian 

if ~isscalar(n) | ~isscalar(l) | ~isnumeric(n) | ~isnumeric(l) 
  error('n, l must be numerical scalars');
end

%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
z1 = [ -9    6*sqrt(3)              0           0           0           0           0;  % 3P
        0    -sqrt(14)     2*sqrt(22)           0           0           0           0;  % 3F
        0            0      sqrt(143)           0           0           0           0;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0]; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
z2 = [ 63   48*sqrt(3)              0           0           0           0           0;  % 3P
        0 112*sqrt(14)    16*sqrt(22)           0           0           0           0;  % 3F
        0            0   -7*sqrt(143)           0           0           0           0;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0]; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
z3 = [  0    1/sqrt(3)              0           0           0           0           0;  % 3P
        0            0    -3/sqrt(22)           0           0           0           0;  % 3F
        0            0              0           0           0           0           0;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0]; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
z4 = [  1            0              0           0           0           0           0;  % 3P
        0            0              0           0           0           0           0;  % 3F
        0            0    9/sqrt(143)           0           0           0           0;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0]; % 1I

z1 = z1+triu(z1,1)'; z2 = z2+triu(z2,1)'; z3 = z3+triu(z3,1)'; z4 = z4+triu(z4,1)';

oldz = {z1 z2 z3 z4}; st_old = racah_states(2,l); s = 1/2; t = 2;

T22.f2 = {z1 z2 z3 z4};

if n==2
  return;
elseif n>(2*(2*l+1))
  error('number of equivalent electrons, n>2(2l+1), greater than allowed orbital angular momenta');
elseif n<2
  error('number of equivalent electrons is invalid (1, 0 or negative!)');
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
    z{im} = zeros(num_states);
    for i = 1:num_states        % Calculates only the upper triangle
      for j = i:num_states      % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
        if abs(S(j)-S(i))<=2 && abs(L(j)-L(i))<=2
          [inz,jnz] = find( cfp{i}'*cfp{j}.*oldz{im} );
          for p = 1:length(inz);
            pi = inz(p); pj = jnz(p);
	    z{im}(i,j) = z{im}(i,j) + ( cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
                                        * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1)) ...
                                        * sixj([S(i) t S(j); Sb(pj) s Sb(pi)]) ...
                                        * sixj([L(i) t L(j); Lb(pj) l Lb(pi)]) * oldz{im}(pi,pj) );
          end
        end % if /\S=0,1,2, /\L=0,1,2
      end
    end
    z{im} = z{im}+z{im}'.*lTriFact;
    z{im} = ((n_calc)/(n_calc-2)) .* z{im}; 
  end

  time = toc; if time<60; display(sprintf('Elapsed time was %0.2g s',time));
                     else display(sprintf('Elapsed time was %0.2g min',time/60)); end;

  %if exist('matname') && ~exist(matfilename,'file')
  %  display(sprintf('Saving %s',matfilename));
  %  save(matfilename,'matrix');
  %end

  T22.(sprintf('%c%1g',lower(racah_lconv(l)),n_calc)) = z;
  
  oldz = z;
  st_old = st_now;

end
