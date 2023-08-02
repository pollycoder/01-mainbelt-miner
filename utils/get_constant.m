%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A file to get all the constants we need. 
% au: Astrononic unit (km)
% mu: Gravity coefficient of sun (km3/s2)
% day:One Earth day (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [au,mu,day]=get_constant()
mu = 1.32712440018e11;                  % Gravity cofficient (km3/s2)
au=149597582.5031;                      % Astronomical unit (km)
day=86400;                              % One Earth day (s)
end
