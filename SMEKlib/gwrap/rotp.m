function y = rotp(p, angle, varargin)
%rotp Rotate point by angle.
%
% (c) 2018 Antti Lehikoinen / Smeklab Ltd

if numel(varargin) == 0
    %rotation around origin
    y = [cos(angle) -sin(angle);sin(angle) cos(angle)]*p;
else
    %rotation around point pc
    pc = varargin{1};
    y = bsxfun(@plus, p, -pc);
    y = [cos(angle) -sin(angle);sin(angle) cos(angle)]*p;
    y = bsxfun(@plus, y, pc);
end

end