%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint function
% No input.
% Object:
%   The integrant
%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=transferEndpoint(input)
q=input.phase.integral;
%q=input.phase.finaltime;
output.objective=q;
end