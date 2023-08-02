%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint function
% No input.
% Object:
%   The integrant
%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=transferEndpoint(input)
intresult=input.phase.integral;
output.object=intresult;
end