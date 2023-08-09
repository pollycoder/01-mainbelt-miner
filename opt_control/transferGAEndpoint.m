%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End point function - all new unit
% Event: GA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=transferGAEndpoint(input)
tf1=input.phase(1).finaltime;
xf1=input.phase(1).finalstate;

% Objective
output.objective=-xf1(7);
end