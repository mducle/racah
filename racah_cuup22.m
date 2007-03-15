function c = racah_cuup22(U,Up)
% Looks up the constant of proportionality c(U,U'(22)), after Racah IV

if ~isvector(U) || length(U)~=2 || ~isvector(Up) || length(Up)~=2
  error('U and Up must be vectors of length 2.');
end

t = [0  0  0  0  0  0  1  0  0;
     0  0  0  0  1  0  0  1  0;
     0  0  1  0  0  1  1  1  0;
     0  0  0  1  1  1  1  1  1;
     0  1  0  1  2  1  0  2  2;
     0  0  1  1  1  2  1  2  1;
     1  0  1  1  0  1  2  1  1;
     0  1  1  1  2  2  1  3  2;
     0  0  0  1  2  1  1  2  2];

      if U==[0 0]; i=1;
  elseif U==[1 0]; i=2;
  elseif U==[1 1]; i=3;
  elseif U==[2 0]; i=4;
  elseif U==[2 1]; i=5;
  elseif U==[3 0]; i=6;
  elseif U==[2 2]; i=7;
  elseif U==[3 1]; i=8;
  elseif U==[4 0]; i=9;
  else; i = 0;
end
      if Up==[0 0]; j=1;
  elseif Up==[1 0]; j=2;
  elseif Up==[1 1]; j=3;
  elseif Up==[2 0]; j=4;
  elseif Up==[2 1]; j=5;
  elseif Up==[3 0]; j=6;
  elseif Up==[2 2]; j=7;
  elseif Up==[3 1]; j=8;
  elseif Up==[4 0]; j=9;
  else; j = 0;
end

if i~=0 && j~=0
  c = t(i,j);
else
  c = 0;
end
