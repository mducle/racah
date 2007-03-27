function [states,st_id] = racah_states(n,l)
% Outputs a list of Racah states with quantum numbers S,L,v,U for configuration l^n

if ischar(n) && ischar(l) && length(n)==2
  state = l;
  l = n(1);
  n = double(n(2))-48;
end

if isnumeric(l)
  switch l
    case 1; l = 'p';
    case 2; l = 'd';
    case 3; l = 'f';
  otherwise
    error('l must be 1,2,3');
  end
end

switch l
  case 'p';
    switch n
      case 1; states = {1/2 'P' 0 [0 0]};                                        % p1	2P
      case 2; states = {{1 'P' 0 [0 0]} {0 'S' 0 [0 0]} {0 'D' 0 [0 0]}};        % p2	3P 1S 1D
      case 3; states = {{3/2 'S' 0 [0 0]} {1/2 'P' 0 [0 0]} {1/2 'D' 0 [0 0]}};  % p3	4S 2P 2D
    otherwise
      error('If l=1, need n=1,2 or 3.');
    end
  case 'd';
    switch n
      case 1; states = { {1/2 'D' 1 [1 0]} };   % d1	2D	1 (10)
      case 2; states = { { 1  'P' 2 [1 1]} ...  % d2	3P	2 (11)
                         { 1  'F' 2 [1 1]} ...  %	3F	2 (11)
                         { 0  'S' 0 [0 0]} ...  %	1S	0 (00)
                         { 0  'D' 2 [2 0]} ...  %	1D	2 (20)
                         { 0  'G' 2 [2 0]} };   %	1G	2 (20)
	      st_id = {'3P' '3F' '1S' '1D' '1G'};
      case 3; states = { {3/2 'P' 3 [1 1]} ...  % d3	4P	3 (11)
                         {3/2 'F' 3 [1 1]} ...  %	4F	3 (11)
                         {1/2 'P' 3 [2 1]} ...  %	2P	3 (21)
                         {1/2 'D' 1 [1 0]} ...  %	2D1	1 (10)
                         {1/2 'D' 3 [2 1]} ...  %	2D2	3 (21)
                         {1/2 'F' 3 [2 1]} ...  %	2F	3 (21)
                         {1/2 'G' 3 [2 1]} ...  %	2G	3 (21)
                         {1/2 'H' 3 [2 1]} };   %	2H	3 (21)
 	      st_id = {'4P' '4F' '2P' '2D1' '2D2' '2F' '2G' '2H'};
      case 4; states = { { 2  'D' 4 [1 0]} ...  % d4	5D	4 (10)
                         { 1  'P' 2 [1 1]} ...  %	3P1	2 (11)
                         { 1  'P' 4 [2 1]} ...  %	3P2	4 (21)
                         { 1  'D' 4 [2 1]} ...  %	3D	4 (21)
                         { 1  'F' 2 [1 1]} ...  %	3F1	2 (11)
                         { 1  'F' 4 [2 1]} ...  %	3F2	4 (21)
                         { 1  'G' 4 [2 1]} ...  %	3G	4 (21)
                         { 1  'H' 4 [2 1]} ...  %	3H	4 (21)
                         { 0  'S' 0 [0 0]} ...  %	1S1	0 (00)
                         { 0  'S' 4 [2 2]} ...  %	1S2	4 (22)
                         { 0  'D' 2 [2 0]} ...  %	1D1	2 (20)
                         { 0  'D' 4 [2 2]} ...  %	1D2	4 (22)
                         { 0  'F' 4 [2 2]} ...  %	1F	4 (22)
                         { 0  'G' 2 [2 0]} ...  %	1G1	2 (20)
                         { 0  'G' 4 [2 2]} ...  %	1G2	4 (22)
                         { 0  'I' 4 [2 2]} };   %	1I	4 (22)
 	      st_id = {'5D' '3P1' '3P2' '3D' '3F1' '3F2' '3G' '3H' '1S1' '1S2' '1D1' '1D2' '1F' '1G1' '1G2' '1I'};
      case 5; states = { {5/2 'S' 5 [0 0]} ...  % d5	6S	5 (00)
                         {3/2 'P' 3 [1 1]} ...  %	4P	3 (11)
                         {3/2 'D' 5 [2 0]} ...  % 	4D	5 (20)
                         {3/2 'F' 3 [1 1]} ...  %	4F	3 (11)
                         {3/2 'G' 5 [2 0]} ...  %	4G	5 (20)
                         {1/2 'S' 5 [2 2]} ...  %	2S	5 (22)
                         {1/2 'P' 3 [2 1]} ...  %	2P	3 (21)
                         {1/2 'D' 1 [1 0]} ...  %	2D1	1 (10)
                         {1/2 'D' 3 [2 1]} ...  %	2D2	3 (21)
                         {1/2 'D' 5 [2 2]} ...  %	2D3	5 (22)
                         {1/2 'F' 3 [2 1]} ...  %	2F1	3 (21)
                         {1/2 'F' 5 [2 2]} ...  %	2F2	5 (22)
                         {1/2 'G' 3 [2 1]} ...  %	2G1	3 (21)
                         {1/2 'G' 5 [2 2]} ...  %	2G2	5 (22)
                         {1/2 'H' 3 [2 1]} ...  %	2H	3 (21)
                         {1/2 'I' 5 [2 2]} };   %	2I	5 (22)
 	      st_id = {'6S' '4P' '4D' '4F' '4G' '2S' '2P' '2D1' '2D2' '2D3' '2F1' '2F2' '2G1' '2G2' '2H' '2I'};
    otherwise
      error('If l=2, need n=1,2,3,4 or 5.');
    end
  case 'f';
    switch n
      case 1; states = { {1/2 'F' 1 [1 0]} };   % f1	2F	1 100 10
      case 2; states = { { 1  'P' 2 [1 1]} ...  % f2	3P	2 110 11
                         { 1  'F' 2 [1 0]} ...  %	3F	2 110 10
                         { 1  'H' 2 [1 1]} ...  %	3H	2 110 11
                         { 0  'S' 0 [0 0]} ...  %	1S	0 000 00
                         { 0  'D' 2 [2 0]} ...  %	1D	2 200 20
                         { 0  'G' 2 [2 0]} ...  %	1G	2 200 20
                         { 0  'I' 2 [2 0]} };   %	1I	2 200 20
 	      st_id = {'2F' '3P' '3F' '3H' '1S' '1D' '1G' '1I'};
      case 3; states = { {3/2 'S' 3 [0 0]} ...  % f3	4S	3 111 00
                         {3/2 'D' 3 [2 0]} ...  %	4D	3 111 20
                         {3/2 'F' 3 [1 0]} ...  %	4F	3 111 10
                         {3/2 'G' 3 [2 0]} ...  %	4G	3 111 20
                         {3/2 'I' 3 [2 0]} ...  %	4I	3 111 20
                         {1/2 'P' 3 [1 1]} ...  %	2P	3 210 11
                         {1/2 'D' 3 [2 0]} ...  %	2D1	3 210 20
                         {1/2 'D' 3 [2 1]} ...  %	2D2	3 210 21
                         {1/2 'F' 1 [1 0]} ...  %	2F1	1 100 10
                         {1/2 'F' 3 [2 1]} ...  %	2F2	3 210 21
                         {1/2 'G' 3 [2 0]} ...  %	2G1	3 210 20
                         {1/2 'G' 3 [2 1]} ...  %	2G2	3 210 21
                         {1/2 'H' 3 [1 1]} ...  %	2H1	3 210 11
                         {1/2 'H' 3 [2 1]} ...  %	2H2	3 210 21
                         {1/2 'I' 3 [2 0]} ...  %	2I	3 210 20
                         {1/2 'K' 3 [2 1]} ...  %	2K	3 210 21
                         {1/2 'L' 3 [2 1]} };   %	2L	3 210 21
 	      st_id = {'4S' '4D' '4F' '4G' '4I' '2P' '2D1' '2D2' '2F1' '2F2' '2G1' '2G2' '2H1' '2H2' '2I' '2K' '2L'};
      case 4; states = { { 2  'S' 4 [0 0]} ...  % f4	5S	4 111 00
                         { 2  'D' 4 [2 0]} ...  %	5D	4 111 20
                         { 2  'F' 4 [1 0]} ...  %	5F	4 111 10
                         { 2  'G' 4 [2 0]} ...  %	5G	4 111 20
                         { 2  'I' 4 [2 0]} ...  %	5I	4 111 20
                         { 1  'P' 2 [1 1]} ...  %	3P1	2 110 11
                         { 1  'P' 4 [1 1]} ...  %	3P2	4 211 11
                         { 1  'P' 4 [3 0]} ...  %	3P3	4 211 30
                         { 1  'D' 4 [2 0]} ...  %	3D1	4 211 20
                         { 1  'D' 4 [2 1]} ...  %	3D2	4 211 21
                         { 1  'F' 2 [1 0]} ...  %	3F1	2 110 10
                         { 1  'F' 4 [1 0]} ...  %	3F2	4 211 10
                         { 1  'F' 4 [2 1]} ...  %	3F3	4 211 21
                         { 1  'F' 4 [3 0]} ...  %	3F4	4 211 30
                         { 1  'G' 4 [2 0]} ...  %	3G1	4 211 20
                         { 1  'G' 4 [2 1]} ...  %	3G2	4 211 21
                         { 1  'G' 4 [3 0]} ...  %	3G3	4 211 30
                         { 1  'H' 2 [1 1]} ...  %	3H1	2 110 11
                         { 1  'H' 4 [1 1]} ...  %	3H2	4 211 11
                         { 1  'H' 4 [2 1]} ...  %	3H3	4 211 21
                         { 1  'H' 4 [3 0]} ...  %	3H4	4 211 30
                         { 1  'I' 4 [2 0]} ...  %	3I1	4 211 20
                         { 1  'I' 4 [3 0]} ...  %	3I2	4 211 30
                         { 1  'K' 4 [2 1]} ...  %	3K1	4 211 21
                         { 1  'K' 4 [3 0]} ...  %	3K2	4 211 30
                         { 1  'L' 4 [2 1]} ...  %	3L	4 211 21
                         { 1  'M' 4 [3 0]} ...  %	3M	4 211 30
                         { 0  'S' 0 [0 0]} ...  %	1S1	0 000 00
                         { 0  'S' 4 [2 2]} ...  %	1S2	4 220 22
                         { 0  'D' 2 [2 0]} ...  %	1D1	2 200 20
                         { 0  'D' 4 [2 0]} ...  %	1D2	4 220 20
                         { 0  'D' 4 [2 1]} ...  %	1D3	4 220 21
                         { 0  'D' 4 [2 2]} ...  %	1D4	4 220 22
                         { 0  'F' 4 [2 1]} ...  %	1F	4 220 21
                         { 0  'G' 2 [2 0]} ...  %	1G1	2 200 20
                         { 0  'G' 4 [2 0]} ...  %	1G2	4 220 20
                         { 0  'G' 4 [2 1]} ...  %	1G3	4 220 21
                         { 0  'G' 4 [2 2]} ...  %	1G4	4 220 22
                         { 0  'H' 4 [2 1]} ...  %	1H1	4 220 21
                         { 0  'H' 4 [2 2]} ...  %	1H2	4 220 22
                         { 0  'I' 2 [2 0]} ...  %	1I1	2 200 20
                         { 0  'I' 4 [2 0]} ...  %	1I2	4 220 20
                         { 0  'I' 4 [2 2]} ...  %	1I3	4 220 22
                         { 0  'K' 4 [2 1]} ...  %	1K	4 220 21
                         { 0  'L' 4 [2 1]} ...  %	1L1	4 220 21
                         { 0  'L' 4 [2 2]} ...  %	1L2	4 220 22
                         { 0  'N' 4 [2 2]} };   %	1N	4 220 22
 	      st_id = {'5S' '5D' '5F' '5G' '5I' '3P1' '3P2' '3P3' '3D1' '3D2' '3F1' '3F2' '3F3' '3F4' '3G1' '3G2' '3G3' '3H1' '3H2' '3H3' '3H4' '3I1' '3I2' '3K1' '3K2' '3L' '3M' '1S1' '1S2' '1D1' '1D2' '1D3' '1D4' '1F' '1G1' '1G2' '1G3' '1G4' '1H1' '1H2' '1I1' '1I2' '1I3' '1K' '1L1' '1L2' '1N'};
      case 5; states = { {5/2 'P' 5 [1 1]} ...  % f5	6P	5 110 11
                         {5/2 'F' 5 [1 0]} ...  %	6F	5 110 10
                         {5/2 'H' 5 [1 1]} ...  %	6H	5 110 11
                         {3/2 'S' 3 [0 0]} ...  %	4S	3 111 00
                         {3/2 'P' 5 [1 1]} ...  %	4P1	5 211 11
                         {3/2 'P' 5 [3 0]} ...  %	4P2	5 211 30
                         {3/2 'D' 3 [2 0]} ...  %	4D1	3 111 20
                         {3/2 'D' 5 [2 0]} ...  %	4D2	5 211 20
                         {3/2 'D' 5 [2 1]} ...  %	4D3	5 211 21
                         {3/2 'F' 3 [1 0]} ...  %	4F1	3 111 10
                         {3/2 'F' 5 [1 0]} ...  %	4F2	5 211 10
                         {3/2 'F' 5 [2 1]} ...  %	4F3	5 211 21
                         {3/2 'F' 5 [3 0]} ...  %	4F4	5 211 30
                         {3/2 'G' 3 [2 0]} ...  %	4G1	3 111 20
                         {3/2 'G' 5 [2 0]} ...  %	4G2	5 211 20
                         {3/2 'G' 5 [2 1]} ...  %	4G3	5 211 21
                         {3/2 'G' 5 [3 0]} ...  %	4G4	5 211 30
                         {3/2 'H' 5 [1 1]} ...  %	4H1	5 211 11
                         {3/2 'H' 5 [2 1]} ...  %	4H2	5 211 21
                         {3/2 'H' 5 [3 0]} ...  %	4H3	5 211 30
                         {3/2 'I' 3 [2 0]} ...  %	4I1	3 111 20
                         {3/2 'I' 5 [2 0]} ...  %	4I2	5 211 20
                         {3/2 'I' 5 [3 0]} ...  %	4I3	5 211 30
                         {3/2 'K' 5 [2 1]} ...  %	4K1	5 211 21
                         {3/2 'K' 5 [3 0]} ...  %	4K2	5 211 30
                         {3/2 'L' 5 [2 1]} ...  %	4L	5 211 21
                         {3/2 'M' 5 [3 0]} ...  %	4M	5 211 30
                         {1/2 'P' 3 [1 1]} ...  %	2P1	3 210 11
                         {1/2 'P' 5 [1 1]} ...  %	2P2	5 211 11
                         {1/2 'P' 5 [3 0]} ...  %	2P3	5 211 30
                         {1/2 'P' 5 [3 1]} ...  %	2P4	5 211 31
                         {1/2 'D' 3 [2 0]} ...  %	2D1	3 210 20
                         {1/2 'D' 3 [2 1]} ...  %	2D2	3 210 21
                         {1/2 'D' 5 [2 0]} ...  %	2D3	5 211 20
                         {1/2 'D' 5 [2 1]} ...  %	2D4 	5 211 21
                         {1/2 'D' 5 [3 1]} ...  %	2D5	5 211 31
                         {1/2 'F' 1 [1 0]} ...  %	2F1	1 100 10
                         {1/2 'F' 3 [2 1]} ...  %	2F2	3 210 21
                         {1/2 'F' 5 [1 0]} ...  %	2F3	5 221 10
                         {1/2 'F' 5 [2 1]} ...  %	2F4	5 221 21
                         {1/2 'F' 5 [3 0]} ...  %	2F5	5 221 30
                         {1/2 'F' 5 [3 1]} ...  %	2F6	5 221 31A
                         {1/2 'Fp' 5 [3 1]} ... %	2F7	5 221 31B
                         {1/2 'G' 3 [2 0]} ...  %	2G1	3 210 20
                         {1/2 'G' 3 [2 1]} ...  %	2G2	3 210 21
                         {1/2 'G' 5 [2 0]} ...  %	2G3	5 221 20
                         {1/2 'G' 5 [2 1]} ...  %	2G4 	5 221 21
                         {1/2 'G' 5 [3 0]} ...  %	2G5	5 221 30
                         {1/2 'G' 5 [3 1]} ...  %	2G6	5 221 31
                         {1/2 'H' 3 [1 1]} ...  %	2H1	3 210 11
                         {1/2 'H' 3 [2 1]} ...  %	2H2	3 210 21
                         {1/2 'H' 5 [1 1]} ...  %	2H3	5 221 11
                         {1/2 'H' 5 [2 1]} ...  %	2H4 	5 221 21
                         {1/2 'H' 5 [3 0]} ...  %	2H5	5 221 30
                         {1/2 'H' 5 [3 1]} ...  %	2H6	5 221 31A
                         {1/2 'Hp' 5 [3 1]} ... %	2H7	5 221 31B
                         {1/2 'I' 3 [2 0]} ...  %	2I1	3 210 20
                         {1/2 'I' 5 [2 0]} ...  %	2I2	5 221 20
                         {1/2 'I' 5 [3 0]} ...  %	2I3	5 221 30
                         {1/2 'I' 5 [3 1]} ...  %	2I4	5 221 31A
                         {1/2 'Ip' 5 [3 1]} ... %	2I5	5 221 31B
                         {1/2 'K' 3 [2 1]} ...  %	2K1	3 210 21
                         {1/2 'K' 5 [2 1]} ...  %	2K2	5 221 21
                         {1/2 'K' 5 [3 0]} ...  %	2K3	5 221 30
                         {1/2 'K' 5 [3 1]} ...  %	2K4	5 221 31A
                         {1/2 'Kp' 5 [3 1]} ... %	2K5	5 221 31B
                         {1/2 'L' 5 [2 1]} ...  %	2L1	3 210 21
                         {1/2 'L' 5 [2 1]} ...  %	2L2	5 221 21
                         {1/2 'L' 5 [3 1]} ...  %	2L3	5 221 31
                         {1/2 'M' 5 [3 0]} ...  %	2M1	5 221 30
                         {1/2 'M' 5 [3 1]} ...  %	2M2	5 221 31
                         {1/2 'N' 5 [3 1]} ...  %	2N	5 221 31
                         {1/2 'O' 5 [3 1]} };   %	2O	5 221 31
  	      st_id = {'6P' '6F' '6H' '4S' '4P1' '4P2' '4D1' '4D2' '4D3' '4F1' '4F2' '4F3' '4F4' '4G1' '4G2' '4G3' '4G4' '4H1' '4H2' '4H3' '4I1' '4I2' '4I3' '4K1' '4K2' '4L' '4M' '2P1' '2P2' '2P3' '2P4' '2D1' '2D2' '2D3' '2D4' '2D5' '2F1' '2F2' '2F3' '2F4' '2F5' '2F6' '2F7' '2G1' '2G2' '2G3' '2G4' '2G5' '2G6' '2H1' '2H2' '2H3' '2H4' '2H5' '2H6' '2H7' '2I1' '2I2' '2I3' '2I4' '2I5' '2K1' '2K2' '2K3' '2K4' '2K5' '2L1' '2L2' '2L3' '2M1' '2M2' '2N' '2O'};
      case 6; states = { { 3  'F' 6 [1 0]} ...  % f6	7F	6 100 10
                         { 2  'S' 4 [0 0]} ...  %	5S	4 111 00
                         { 2  'P' 6 [1 1]} ...  %	5P	6 210 11
                         { 2  'D' 4 [2 0]} ...  %	5D1	4 111 20
                         { 2  'D' 6 [2 0]} ...  %	5D2	6 210 20
                         { 2  'D' 6 [2 1]} ...  %	5D3	6 210 21
                         { 2  'F' 4 [1 0]} ...  %	5F1	4 111 10
                         { 2  'F' 6 [2 1]} ...  %	5F2	6 210 21
                         { 2  'G' 4 [2 0]} ...  %	5G1	4 111 20
                         { 2  'G' 6 [2 0]} ...  %	5G2	6 210 20
                         { 2  'G' 6 [2 1]} ...  %	5G3	6 210 21
                         { 2  'H' 6 [1 1]} ...  %	5H1	6 210 11
                         { 2  'H' 6 [2 1]} ...  %	5H2	6 210 21
                         { 2  'I' 4 [2 0]} ...  %	5I1	4 111 20
                         { 2  'I' 6 [2 0]} ...  %	5I2	6 210 20
                         { 2  'K' 6 [2 1]} ...  %	5K	6 210 21
                         { 2  'L' 6 [2 1]} ...  %	5L	6 210 21
                         { 1  'P' 2 [1 1]} ...  %	3P1	2 110 11
                         { 1  'P' 4 [1 1]} ...  %	3P2	4 211 11
                         { 1  'P' 4 [3 0]} ...  %	3P3	4 211 30
                         { 1  'P' 6 [1 1]} ...  %	3P4	6 221 11
                         { 1  'P' 6 [3 0]} ...  %	3P5	6 221 30
                         { 1  'P' 6 [3 1]} ...  %	3P6	6 221 31
                         { 1  'D' 4 [2 0]} ...  %	3D1	4 211 20
                         { 1  'D' 4 [2 1]} ...  %	3D2	4 211 21
                         { 1  'D' 6 [2 0]} ...  %	3D3	6 221 20
                         { 1  'D' 6 [2 1]} ...  %	3D4	6 221 21
                         { 1  'D' 6 [3 1]} ...  %	3D5	6 221 31
                         { 1  'F' 2 [1 0]} ...  %	3F1	2 110 10
                         { 1  'F' 4 [1 0]} ...  %	3F1	4 211 10
                         { 1  'F' 4 [2 1]} ...  %	3F3	4 211 21
                         { 1  'F' 4 [3 0]} ...  %	3F4	4 211 30
                         { 1  'F' 6 [1 0]} ...  %	3F5	6 221 10
                         { 1  'F' 6 [2 1]} ...  %	3F6	6 221 21
                         { 1  'F' 6 [3 0]} ...  %	3F7	6 221 30
                         { 1  'F' 6 [3 1]} ...  %	3F8	6 221 31A
                         { 1  'Fp' 6 [3 1]} ... %	3F9	6 221 31B
                         { 1  'G' 4 [2 0]} ...  %	3G1	4 211 20
                         { 1  'G' 4 [2 1]} ...  %	3G2	4 211 21
                         { 1  'G' 4 [3 0]} ...  %	3G3	4 211 30
                         { 1  'G' 6 [2 0]} ...  %	3G4	6 221 20
                         { 1  'G' 6 [2 1]} ...  %	3G5	6 221 21
                         { 1  'G' 6 [3 0]} ...  %	3G6	6 221 30
                         { 1  'G' 6 [3 1]} ...  %	3G7	6 221 31
                         { 1  'H' 2 [1 1]} ...  %	3H1	2 110 11
                         { 1  'H' 4 [1 1]} ...  %	3H2	4 211 11
                         { 1  'H' 4 [2 1]} ...  %	3H3	4 211 21
                         { 1  'H' 4 [3 0]} ...  %	3H4	4 211 30
                         { 1  'H' 6 [1 1]} ...  %	3H5	6 221 11
                         { 1  'H' 6 [2 1]} ...  %	3H6	6 221 21
                         { 1  'H' 6 [3 0]} ...  %	3H7	6 221 30
                         { 1  'H' 6 [3 1]} ...  %	3H8	6 221 31A
                         { 1  'Hp' 6 [3 1]} ... %	3H9	6 221 31B
                         { 1  'I' 4 [2 0]} ...  %	3I1	4 211 20
                         { 1  'I' 4 [3 0]} ...  %	3I2	4 211 30
                         { 1  'I' 6 [2 0]} ...  %	3I3	6 221 20
                         { 1  'I' 6 [3 0]} ...  %	3I4	6 221 30
                         { 1  'I' 6 [3 1]} ...  %	3I5	6 221 31A
                         { 1  'Ip' 6 [3 1]} ... %	3I6	6 221 31B
                         { 1  'K' 4 [2 1]} ...  %	3K1	4 211 21
                         { 1  'K' 4 [3 0]} ...  %	3K2	4 211 30
                         { 1  'K' 6 [2 1]} ...  %	3K3	6 221 21
                         { 1  'K' 6 [3 0]} ...  %	3K4	6 221 30
                         { 1  'K' 6 [3 1]} ...  %	3K5	6 221 31A
                         { 1  'Kp' 6 [3 1]} ... %	3K6	6 221 31B
                         { 1  'L' 4 [2 1]} ...  %	3L1	4 211 21
                         { 1  'L' 6 [2 1]} ...  %	3L2	6 221 21
                         { 1  'L' 6 [3 1]} ...  %	3L3	6 221 31
                         { 1  'M' 4 [3 0]} ...  %	3M1	4 211 30
                         { 1  'M' 6 [3 0]} ...  %	3M2	6 221 30
                         { 1  'M' 6 [3 1]} ...  %	3M3	6 221 31
                         { 1  'N' 6 [3 1]} ...  %	3N	6 221 31
                         { 1  'O' 6 [3 1]} ...  %	3O	6 221 31
                         { 0  'S' 0 [0 0]} ...  %	1S1	0 000 00
                         { 0  'S' 4 [2 2]} ...  %	1S2	4 220 22
                         { 0  'S' 6 [0 0]} ...  %	1S3	6 222 00
                         { 0  'S' 6 [4 0]} ...  %	1S4	6 222 40
                         { 0  'P' 6 [3 0]} ...  %	1P	6 222 30
                         { 0  'D' 2 [2 0]} ...  %	1D1	2 200 20
                         { 0  'D' 4 [2 0]} ...  %	1D2	4 220 20
                         { 0  'D' 4 [2 1]} ...  %	1D3	4 220 21
                         { 0  'D' 4 [2 2]} ...  %	1D4	4 220 22
                         { 0  'D' 6 [2 0]} ...  %	1D5	6 222 20
                         { 0  'D' 6 [4 0]} ...  %	1D6	6 222 40
                         { 0  'F' 4 [2 1]} ...  %	1F1	4 220 21
                         { 0  'F' 6 [1 0]} ...  %	1F2	6 222 10
                         { 0  'F' 6 [3 0]} ...  %	1F3	6 222 30
                         { 0  'F' 6 [4 0]} ...  %	1F4	6 222 40
                         { 0  'G' 2 [2 0]} ...  %	1G1	2 200 20
                         { 0  'G' 4 [2 0]} ...  %	1G2	4 220 20
                         { 0  'G' 4 [2 1]} ...  %	1G3	4 220 21
                         { 0  'G' 4 [2 2]} ...  %	1G4	4 220 22
                         { 0  'G' 6 [2 0]} ...  %	1G5	6 222 20
                         { 0  'G' 6 [3 0]} ...  %	1G6	6 222 30
                         { 0  'G' 6 [4 0]} ...  %	1G7	6 222 40A
                         { 0  'G' 6 [4 0]} ...  %	1G8	6 222 40B
                         { 0  'H' 4 [2 1]} ...  %	1H1	4 220 21
                         { 0  'H' 4 [2 2]} ...  %	1H2	4 220 22
                         { 0  'H' 6 [3 0]} ...  %	1H3	6 222 30
                         { 0  'H' 6 [4 0]} ...  %	1H4	6 222 40
                         { 0  'I' 2 [2 0]} ...  %	1I1	2 200 20
                         { 0  'I' 4 [2 0]} ...  %	1I2	4 220 20
                         { 0  'I' 4 [2 2]} ...  %	1I3	4 220 22
                         { 0  'I' 6 [2 0]} ...  %	1I4	6 222 20
                         { 0  'I' 6 [3 0]} ...  %	1I5	6 222 30
                         { 0  'I' 6 [4 0]} ...  %	1I6	6 222 40A
                         { 0  'Ip' 6 [4 0]} ... %	1I7	6 222 40B
                         { 0  'K' 4 [2 1]} ...  %	1K1	4 220 21
                         { 0  'K' 6 [3 0]} ...  %	1K2	6 222 30
                         { 0  'K' 6 [4 0]} ...  %	1K3	6 222 40
                         { 0  'L' 4 [2 1]} ...  %	1L1	4 220 21
                         { 0  'L' 4 [2 2]} ...  %	1L2	4 220 22
                         { 0  'L' 6 [4 0]} ...  %	1L3	6 222 40A
                         { 0  'Lp' 6 [4 0]} ... %	1L4	6 222 40B
                         { 0  'M' 6 [3 0]} ...  %	1M1	6 222 30
                         { 0  'M' 6 [4 0]} ...  %	1M2	6 222 40
                         { 0  'N' 4 [2 2]} ...  %	1N1	4 220 22
                         { 0  'N' 6 [4 0]} ...  %	1N2	6 222 40
                         { 0  'Q' 6 [4 0]} };   %	1Q	6 222 40
	      st_id = {'7F' '5S' '5P' '5D1' '5D2' '5D3' '5F1' '5F2' '5G1' '5G2' '5G3' '5H1' '5H2' '5I1' '5I2' '5K' '5L' '3P1' '3P2' '3P3' '3P4' '3P5' '3P6' '3D1' '3D2' '3D3' '3D4' '3D5' '3F1' '3F1' '3F3' '3F4' '3F5' '3F6' '3F7' '3F8' '3F9' '3G1' '3G2' '3G3' '3G4' '3G5' '3G6' '3G7' '3H1' '3H2' '3H3' '3H4' '3H5' '3H6' '3H7' '3H8' '3H9' '3I1' '3I2' '3I3' '3I4' '3I5' '3I6' '3K1' '3K2' '3K3' '3K4' '3K5' '3K6' '3L1' '3L2' '3L3' '3M1' '3M2' '3M3' '3N' '3O' '1S1' '1S2' '1S3' '1S4' '1P' '1D1' '1D2' '1D3' '1D4' '1D5' '1D6' '1F1' '1F2' '1F3' '1F4' '1G1' '1G2' '1G3' '1G4' '1G5' '1G6' '1G7' '1G8' '1H1' '1H2' '1H3' '1H4' '1I1' '1I2' '1I3' '1I4' '1I5' '1I6' '1I7' '1K1' '1K2' '1K3' '1L1' '1L2' '1L3' '1L4' '1M1' '1M2' '1N1' '1N2' '1Q' '8S'};
      case 7; states = { {7/2 'S' 7 [0 0]} ...  % 7f	8S	7 000 00
                         {5/2 'P' 5 [1 1]} ...  %	6P	5 110 11
                         {5/2 'D' 7 [2 0]} ...  %	6D	7 200 20
                         {5/2 'F' 5 [1 0]} ...  %	6F	5 110 10
                         {5/2 'G' 7 [2 0]} ...  %	6G	7 200 20
                         {5/2 'H' 5 [1 1]} ...  %	6H	5 110 11
                         {5/2 'I' 7 [2 0]} ...  %	6I	7 200 20
                         {3/2 'S' 3 [0 0]} ...  %	4S1	3 111 00
                         {3/2 'S' 7 [2 2]} ...  %	4S2	7 220 22
                         {3/2 'P' 5 [1 1]} ...  %	4P1	5 211 11
                         {3/2 'P' 5 [3 0]} ...  %	4P2	5 211 30
                         {3/2 'D' 3 [2 0]} ...  %	4D1	3 111 20
                         {3/2 'D' 5 [2 0]} ...  %	4D2	5 211 20
                         {3/2 'D' 5 [2 1]} ...  %	4D3	5 211 21
                         {3/2 'D' 7 [2 0]} ...  %	4D4	7 220 20
                         {3/2 'D' 7 [2 1]} ...  %	4D5	7 220 21
                         {3/2 'D' 7 [2 2]} ...  %	4D6	7 220 22
                         {3/2 'F' 3 [1 0]} ...  %	4F1	3 111 10
                         {3/2 'F' 5 [1 0]} ...  %	4F2	5 211 10
                         {3/2 'F' 5 [2 1]} ...  %	4F3	5 211 21
                         {3/2 'F' 5 [3 0]} ...  %	4F4	5 211 30
                         {3/2 'F' 7 [2 1]} ...  %	4F5	7 220 21
                         {3/2 'G' 3 [2 0]} ...  %	4G1	3 111 20
                         {3/2 'G' 5 [2 0]} ...  %	4G2	5 211 20
                         {3/2 'G' 5 [2 1]} ...  %	4G3	5 211 21
                         {3/2 'G' 5 [3 0]} ...  %	4G4	5 211 30
                         {3/2 'G' 7 [2 0]} ...  %	4G5	7 220 20
                         {3/2 'G' 7 [2 1]} ...  %	4G6	7 220 21
                         {3/2 'G' 7 [2 2]} ...  %	4G7	7 220 22
                         {3/2 'H' 5 [1 1]} ...  %	4H1	5 211 11
                         {3/2 'H' 5 [2 1]} ...  %	4H2	5 211 21
                         {3/2 'H' 5 [3 0]} ...  %	4H3	5 211 30
                         {3/2 'H' 7 [2 1]} ...  %	4H4	7 220 21
                         {3/2 'H' 7 [2 2]} ...  %	4H5	7 220 22
                         {3/2 'I' 3 [2 0]} ...  %	4I1	3 111 20
                         {3/2 'I' 5 [2 0]} ...  %	4I2	5 211 20
                         {3/2 'I' 5 [3 0]} ...  %	4I3	5 211 30
                         {3/2 'I' 7 [2 0]} ...  %	4I4	7 220 20
                         {3/2 'I' 7 [2 2]} ...  %	4I5	7 220 22
                         {3/2 'K' 5 [2 1]} ...  %	4K1	5 211 21
                         {3/2 'K' 5 [3 0]} ...  %	4K2	5 211 30
                         {3/2 'K' 7 [2 1]} ...  %	4K3	7 220 21
                         {3/2 'L' 5 [2 1]} ...  %	4L1	5 211 21
                         {3/2 'L' 7 [2 1]} ...  %	4L2	7 220 21
                         {3/2 'L' 7 [2 2]} ...  %	4L3	7 220 22
                         {3/2 'M' 5 [3 0]} ...  %	4M	5 211 30
                         {3/2 'N' 7 [2 2]} ...  %	4N	7 220 22
                         {1/2 'S' 7 [0 0]} ...  %	2S1	7 222 00
                         {1/2 'S' 7 [4 0]} ...  %	2S2	7 222 40
                         {1/2 'P' 3 [1 1]} ...  %	2P1	3 210 11
                         {1/2 'P' 5 [1 1]} ...  %	2P2	5 221 11
                         {1/2 'P' 5 [3 0]} ...  %	2P3	5 221 30
                         {1/2 'P' 5 [3 1]} ...  %	2P4	5 221 31
                         {1/2 'P' 7 [3 0]} ...  %	2P5	7 222 30
                         {1/2 'D' 3 [2 0]} ...  %	2D1	3 210 20
                         {1/2 'D' 3 [2 1]} ...  %	2D2	3 210 21
                         {1/2 'D' 5 [2 0]} ...  %	2D3	5 221 20
                         {1/2 'D' 5 [2 1]} ...  %	2D4	5 221 21
                         {1/2 'D' 5 [3 1]} ...  %	2D5	5 221 31
                         {1/2 'D' 7 [2 0]} ...  %	2D6	7 222 20
                         {1/2 'D' 7 [4 0]} ...  %	2D7	7 222 40
                         {1/2 'F' 1 [1 0]} ...  %	2F1	1 100 10
                         {1/2 'F' 3 [2 1]} ...  %	2F2	3 210 21
                         {1/2 'F' 5 [1 0]} ...  %	2F3	5 221 10
                         {1/2 'F' 5 [2 1]} ...  %	2F4	5 221 21
                         {1/2 'F' 5 [3 0]} ...  %	2F5	5 221 30
                         {1/2 'F' 5 [3 1]} ...  %	2F6	5 221 31A
                         {1/2 'Fp' 5 [3 1]} ... %	2F7	5 221 31B
                         {1/2 'F' 7 [1 0]} ...  %	2F8	7 222 10
                         {1/2 'F' 7 [3 0]} ...  %	2F9	7 222 30
                         {1/2 'F' 7 [4 0]} ...  %	2F10	7 222 40
                         {1/2 'G' 3 [2 0]} ...  %	2G1	3 210 20
                         {1/2 'G' 3 [2 1]} ...  %	2G2	3 210 21
                         {1/2 'G' 5 [2 0]} ...  %	2G3	5 221 20
                         {1/2 'G' 5 [2 1]} ...  %	2G4	5 221 21
                         {1/2 'G' 5 [3 0]} ...  %	2G5	5 221 30
                         {1/2 'G' 5 [3 1]} ...  %	2G6	5 221 31
                         {1/2 'G' 7 [2 0]} ...  %	2G7	7 222 20
                         {1/2 'G' 7 [3 0]} ...  %	2G8	7 222 30
                         {1/2 'G' 7 [4 0]} ...  %	2G9	7 222 40A
                         {1/2 'Gp' 7 [4 0]} ... %	2G10	7 222 40B
                         {1/2 'H' 3 [1 1]} ...  %	2H1	3 210 11
                         {1/2 'H' 3 [2 1]} ...  %	2H2	3 210 21
                         {1/2 'H' 5 [1 1]} ...  %	2H3	5 221 11
                         {1/2 'H' 5 [2 1]} ...  %	2H4	5 221 21
                         {1/2 'H' 5 [3 0]} ...  %	2H5	5 221 30
                         {1/2 'H' 5 [3 1]} ...  %	2H6	5 221 31A
                         {1/2 'Hp' 5 [3 1]} ... %	2H7	7 222 31B
                         {1/2 'H' 7 [3 0]} ...  %	2H8	7 222 30
                         {1/2 'H' 7 [4 0]} ...  %	2H9	7 222 40
                         {1/2 'I' 3 [2 0]} ...  %	2I1	3 210 20
                         {1/2 'I' 5 [2 0]} ...  %	2I2	5 221 20
                         {1/2 'I' 5 [3 0]} ...  %	2I3	5 221 30
                         {1/2 'I' 5 [3 1]} ...  %	2I4	5 221 31A
                         {1/2 'Ip' 5 [3 1]} ... %	2I5	5 221 31B
                         {1/2 'I' 7 [2 0]} ...  %	2I6	7 222 20
                         {1/2 'I' 7 [3 0]} ...  %	2I7	7 222 30
                         {1/2 'I' 7 [4 0]} ...  %	2I8	7 222 40A
                         {1/2 'Ip' 7 [4 0]} ... %	2I9	7 222 40B
                         {1/2 'K' 3 [2 1]} ...  %	2K1	3 210 21
                         {1/2 'K' 5 [2 1]} ...  %	2K2	5 221 21
                         {1/2 'K' 5 [3 0]} ...  %	2K3	5 221 30
                         {1/2 'K' 5 [3 1]} ...  %	2K4	5 221 31A
                         {1/2 'Kp' 5 [3 1]} ... %	2K5	5 221 31B
                         {1/2 'K' 7 [3 0]} ...  %	2K6	7 222 30
                         {1/2 'K' 7 [4 0]} ...  %	2K7	7 222 40
                         {1/2 'L' 3 [2 1]} ...  %	2L1	3 210 21
                         {1/2 'L' 5 [2 1]} ...  %	2L2	5 221 21
                         {1/2 'L' 5 [3 1]} ...  %	2L3	5 221 31
                         {1/2 'L' 7 [4 0]} ...  %	2L4	7 222 40A
                         {1/2 'Lp' 7 [4 0]} ... %	2L5	7 222 40B
                         {1/2 'M' 5 [3 0]} ...  %	2M1	5 221 30
                         {1/2 'M' 5 [3 1]} ...  %	2M2	5 221 31
                         {1/2 'M' 7 [3 0]} ...  %	2M3	7 222 30
                         {1/2 'M' 7 [4 0]} ...  %	2M4	7 222 40
                         {1/2 'N' 5 [3 1]} ...  %	2N1	5 221 31
                         {1/2 'N' 7 [4 0]} ...  %	2N2	7 222 40
                         {1/2 'O' 5 [3 1]} ...  %	2O	5 221 31
                         {1/2 'Q' 7 [4 0]} };   %	2Q	7 222 40
 	      st_id = {'8S' '6P' '6D' '6F' '6G' '6H' '6I' '4S1' '4S2' '4P1' '4P2' '4D1' '4D2' '4D3' '4D4' '4D5' '4D6' '4F1' '4F2' '4F3' '4F4' '4F5' '4G1' '4G2' '4G3' '4G4' '4G5' '4G6' '4G7' '4H1' '4H2' '4H3' '4H4' '4H5' '4I1' '4I2' '4I3' '4I4' '4I5' '4K1' '4K2' '4K3' '4L1' '4L2' '4L3' '4M' '4N' '2S1' '2S2' '2P1' '2P2' '2P3' '2P4' '2P5' '2D1' '2D2' '2D3' '2D4' '2D5' '2D6' '2D7' '2F1' '2F2' '2F3' '2F4' '2F5' '2F6' '2F7' '2F8' '2F9' '2F10' '2G1' '2G2' '2G3' '2G4' '2G5' '2G6' '2G7' '2G8' '2G9' '2G10' '2H1' '2H2' '2H3' '2H4' '2H5' '2H6' '2H7' '2H8' '2H9' '2I1' '2I2' '2I3' '2I4' '2I5' '2I6' '2I7' '2I8' '2I9' '2K1' '2K2' '2K3' '2K4' '2K5' '2K6' '2K7' '2L1' '2L2' '2L3' '2L4' '2L5' '2M1' '2M2' '2M3' '2M4' '2N1' '2N2' '2O' '2Q'};
    otherwise
      error('If l=3, need n=1,2,3,4,5,6 or 7.');
    end
  otherwise
    error('l must be p,d,f');
end

if exist('state')
  for i = 1:length(states)
    if strcmp(st_id{i},state)
      break;
    end
  end
  states = states{i};
end
