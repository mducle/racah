function emat = racah_emat(n,l,E,states)
% Calculates the electrostatic matrix after the method of Racah. 
%
% Syntax:  emat = racah_emat(n,l,E)
%
% Inputs:  n    - scalar - the number of electrons in the outer shell
%          l    - scalar - the angular momentum of a single electron in the outer shell
%                          at present only the l=3 (f-electron) shell is implemented.
%          E    - vector - [E0 E1 E2 E3] - the Racah electrostatic parameters, a linear
%                          combination of the Slater Integrals F1,F2,F3. Use the function
%                          racah_FtoE and racah_EtoF to convert between the two parameters.
% Outputs: emat - matrix - the energy matrix in terms of the states of the l^n 
%                          configurations. These are listed by the function racah_states.

% References: G. Racah, "Theory of Complex Spectra IV", Phys. Rev., vol. 76, pp1352 (1949)

% By Duc Le - 2007 - duc.le@ucl.ac.uk

% Sun Feb 25 23:43:13 GMT 2007 - initial release.
% mdl - 070303 - updated to take arbitrary states, rather than all states of a f^n conf.

if ~exist('states')
  states = racah_states(n,l);
end

e0   = n*(n-1)/2;

for i = 1:length(states)
  v = states{i}{3}; S = states{i}{1};
  e1(i,i) = (9/2)*(n-v) + v*(v+2)/4 - S*(S+1);
end

emat = (E(1)*e0).*eye(length(states)) ...
      + E(2).*e1 ...
      + E(3).*racah_e2(n,l,states) ...
      + E(4).*racah_e3(n,l,states);
