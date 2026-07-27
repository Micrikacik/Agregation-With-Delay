function [intRad] = getIntRad(d)

base2DIntRad = 0.05;                                    % interaction radius in 2D
kappa_2 = pi;                                           % "volume" of a unit 2-ball
kappa = pi^(d/2) / gamma(d / 2 + 1);                    % volume of a unit d-ball
intRad = (kappa_2 / kappa * base2DIntRad^2)^(1/d);      % interaction radius in d-dim space