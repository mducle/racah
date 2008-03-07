function [A2,A4,A6] = wy2stev(B2,B4,B6,no)
% Converts the crystal field parameters in Wybourne normalisation to Steven's normalisation
%
% Syntax:  [A2,A4,A6] = wy2stev(B2,B4,B6)
%
% Inputs:  B2 - vector - a 1x5 vector of the rank 2 parameters
%             - cell   - a {[1x5] [1x9] [1x13]} cell array of the CF parameters
%          B4 - vector - a 1x9 vector of the rank 4 parameters
%             - string - either 'w' or 's' denoting the normalisation of the input parameters
%                        being either Wybourne or Stevens respectively.
%          B6 - vector - a 1x13 vector of the rank 6 parameters
%          no - string - either 'w' or 's' denoting the normalisation of the input parameters
%
% Outputs: A2 - vector - a 1x5 vector of the converted rank 2 parameters
%             - cell   - a {[1x5] [1x9] [1x13]} cell array of the converted CF parameters
%          A4 - vector - a 1x9 vector of the converted rank 4 parameters
%          A6 - vector - a 1x13 vector of the converted rank 6 parameters

if (exist('B4') && isstr(B4) && strncmp(B4,'s',1)) | (exist('no') && ~strncmp(no,'w',1))
  stevflag = 1;     
end

if iscell(B2) && length(B2)==3
  B4 = B2{2};
  B6 = B2{3};
  B2 = B2{1};
end

conv2 = [2/sqrt(-6)    -1/sqrt(-6)  2  -1/sqrt(6)  2/sqrt(6)];
conv4 = [4/sqrt(-10)   -2/sqrt(-5)  8  -2/sqrt(5)  4/sqrt(10)];
conv6 = [16/sqrt(-105) -8/sqrt(-42) 16 -8/sqrt(42) 16/sqrt(105)];

conv4 = [8/sqrt(-70)    -2/sqrt(-35)  conv4 -2/sqrt(35)  8/sqrt(70)];
conv6 = [16/3/sqrt(-14) -8/sqrt(-105) conv6 -8/sqrt(105) 16/3/sqrt(14)];

conv6 = [16/sqrt(-231) -8/3/sqrt(-77) conv6 -8/3/sqrt(77) 16/sqrt(231)];

if nargout==3
  if exist('stevflag')
    A2 = B2.*conv2;
    A4 = B4.*conv4;
    A6 = B6.*conv6;
  else
    A2 = B2./conv2;
    A4 = B4./conv4;
    A6 = B6./conv6;
  end
else
  if exist('stevflag')
    A2{1} = B2.*conv2;
    A2{2} = B4.*conv4;
    A2{3} = B6.*conv6;
  else
    A2{1} = B2./conv2;
    A2{2} = B4./conv4;
    A2{3} = B6./conv6;
  end
end

