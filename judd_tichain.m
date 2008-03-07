function ti_out = judd_tichain(n,l,T)
% Calculates and saves all the matrices ti for i=2,3,4,6,7,8 (after Judd 1966) up to f^7 
% by a chain calculation and saves it in a file 'timat.mat'

% By Duc Le 2007 - duc.le@ucl.ac.uk

if exist('T') && isnumeric(T) 
  if length(T)==6
    T = [0 T(1:3) 0 T(4:6) 0 0];
  elseif length(T)>8
  else
    error('T must be a vector of length 6, 8 or 10');
  end
end

if exist('timat.mat','file')==2
  load 'timat.mat';
else
  ti.f3 = judd_tif3(); 
end

if ~isnumeric(n) || ~isscalar(n) || n>(2*(2*l+1))
  error('n must be a numeric scalar less than 2(2l+1)');
elseif n==3
  ti_out = ti.f3;
  if exist('T'); ti_out = sparse(17,17); for i = [2 3 4 6 7 8]; ti_out = ti_out + T(i).*ti.f3{i}; end; end
  return;
elseif n<3
  error('n must be >= 3');
end

confstring = sprintf('f%1g',n);
if ~isfield(ti,confstring)
  tic; display('Calculating t_i');
  for i = [2 3 4 6 7 8]
    time_s = toc;
    for in = 4:n
      time_in = toc;
      ti.(sprintf('f%1g',in)){i} = racah_chaincal(in,3,in-1,ti.(sprintf('f%1g',in-1)){i},0);
      if (toc-time_in)>60; display('Calculated t_%1g for f^%1g in %0.2g min',i,n,(toc-time_in)/60); end
    end
    if (toc-time_s)>60; display(sprintf('Calculated t_%1g in %0.2g min',i,toc/60)); end
  end
  if toc>60; display(sprintf('Time Elapsed: %0.2g min\n',i,toc/60));
        else display(sprintf('Time Elapsed: %0.2g s\n',i,toc)); end
end

if exist('T')
  ti_out = zeros(size(ti.(confstring){2}));
  for i = [2 3 4 6 7 8]
    ti_out = ti_out + T(i) .* ti.(confstring){i};
  end
  ti_out = sparse(ti_out);
else
  ti_out = ti.(confstring);
end

save('timat.mat','ti');

%for i = [2 3 4 6 7 8]
%  display(sprintf('Calculating t_%1g',i));
%  for n = 4:7
%    varname = sprintf('t%1g_f%1g',i,n);
%    s.(varname) = racah_chaincal(n,3,n-1,ti{i});
%    ti{i} = s.(varname);
%  end
%  display(sprintf('\n'));
%end
%
%save('timat.mat','-struct','s');
%save('timats.mat','s');

% Total elapsed time approx 5min on a Pentium M 1.6GHz with 756MB memory
