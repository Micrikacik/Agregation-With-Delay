function [intRad] = getIntRad(d)

% Calculates the radius of a d-ball, such that its d-volume is equal to the
% area of a circle with a radius 0.05
% 
% INPUT:
%   d (positive integer array) - dimension/s, for which to calculate the
%       radius
%
% OUTPUT:
%   intRad (positive float array) - radius of a d-ball with the volume 
%       equal to the area of a circle with a radius 0.05, its dimensions
%       are the same as of the input 'd'

arguments
    d double {mustBeInteger, mustBePositive}
end

base2DIntRad = 0.05;                                    % interaction radius in 2D
kappa_2 = pi;                                           % "volume" of a unit 2-ball
kappa = pi.^(d/2) ./ gamma(d / 2 + 1);                    % volume of a unit d-ball
intRad = (kappa_2 ./ kappa * base2DIntRad^2).^(1./d);      % interaction radius in d-dim space