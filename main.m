clc
clear all
%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% Two-body problem
%%%%%%%%%%%%%%%%%%%
range=[0,5;0,5];
[X,result]=pso_parallel(@(x)obj_func(x),range,10000,1000);


