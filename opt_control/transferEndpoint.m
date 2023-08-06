%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint function - all new unit
% No input.
% Object:
%   The integrant
%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=transferEndpoint(input)
%q=input.phase.integral;
%q=input.phase.finaltime;
q=input.phase.finalstate(:,7);
%p=input.phase.initialstate(:,7);
output.objective=-q;
end