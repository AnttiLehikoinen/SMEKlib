function d = angleDifference(a1, a2)
%angleDifference signed difference between angles
% 
% d = angleDifference(a1, a2)
% calculates the signed difference a1-a2 between the angles
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

d = mod(a1 - a2 + pi,2*pi) - pi;
end