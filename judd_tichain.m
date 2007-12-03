% Calculates and saves all the matrices ti for i=2,3,4,6,7,8 (after Judd 1966) up to f^7 
% by a chain calculation and saves it in a file 'timat.mat'

% By Duc Le 2007 - duc.le@ucl.ac.uk

ti = judd_tif3(); 

for i = [2 3 4 6 7 8]
  display(sprintf('Calculating t_%1g',i));
  for n = 4:7
    varname = sprintf('t%1g_f%1g',i,n);
    s.(varname) = racah_chaincal(n,3,n-1,ti{i});
    ti{i} = s.(varname);
  end
  display(sprintf('\n'));
end

save('timat.mat','-struct','s');
save('timats.mat','s');

% Total elapsed time approx 5min on a Pentium M 1.6GHz with 756MB memory
