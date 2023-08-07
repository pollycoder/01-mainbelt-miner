%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A file to get all the constants we need. 
% au: Astrononic unit (km)
% mu: Gravity coefficient of sun (km3/s2)
% day:One Earth day (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [au,mu,day,year]=get_constant()
au=149597582503.1;                      % Astronomical unit (m)
day=86400;                              % One Earth day (s)
mu=1.32712440018e20;                    % Gravity cofficient (m3/s2)
year=2*pi*sqrt(au^3/mu);                % Earth year (s)
end
