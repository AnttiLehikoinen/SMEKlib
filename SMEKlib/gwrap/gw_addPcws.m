function this = gw_addPcws(this, varargin)
%gw_addPcws add surface defined piecewise by segments
% 
% Call syntax
% gw_addPcws(gw, segment1, segments2, ..., surfaces_inside) or
% gw_addPcws(gw, segment1, segments2, ..., isHole)
%
% Each segment is defined by one of the following sets of arguments:
%
% line: 'line', xstart, xend, tol
%
% arc: 'arc', center, astart, aend, tol
%   where astart and aend can be either scalars corresponding to the polar
%   coordinate, or 2x1 vectors describing a point in the xy-plane
%
% arc: 'arc', center, astart, aend, radius, tol
%
% Each line segment consists of >=2 nodes, with mutual distances <= tol.
%
% (c) 2017 Antti Lehikoinen / Aalto University

%array for holding the surface points
x = zeros(3, 1000); np = 0;
function addx(xp)
    % nested function for incrementing the array
    if (np+size(xp,2)) > size(x,2)
        nadd = min(2*size(x,2), size(xp,2));
        x = [x zeros(3,nadd)];
    end
    x(1:size(xp,1),(np+1):(np+size(xp,2))) = xp;
    np = np + size(xp,2);
end

%array for holding named lines, if any
lnames = cell(1, 1000); nln = 0;
function addl(lines, name)
    if (nln + 2) < size(lnames, 2)
        lnames = [lnames cell(1,1000)];
    end
    lnames{nln+1} = np - numel(lines) + lines;
    lnames{nln+2} = name;
    nln = nln + 2;
end

%going through the input until end is reached
ri = 1;
while true
    if strcmpi(varargin{ri}, 'line')
        % line segment
        xstart = varargin{ri+1};
        xend = varargin{ri+2};
      
        tol = varargin{ri+3};
        ri = ri + 4;
        
        Np = ceil( norm(xend(1:2)-xstart(1:2))/tol ) + 1;
        xp = bsxfun(@plus, xstart(1:2,:), bsxfun(@times, xend(1:2,:)-xstart(1:2,:), linspace(0,1,Np)));
        
        if numel(xstart)==3
            xp(3,1) = xstart(3);
        end
        if numel(xend)==3
            xp(3,end) = xend(3);
        end
        
        addx(xp(:,1:(end-1)));
        
        if strcmpi(varargin{ri}, 'linename')
            addl(1:(size(xp,2)-1), varargin{ri+1});
            ri = ri + 2;
        end        
    elseif strcmpi(varargin{ri}, 'arc')
        center = varargin{ri+1};
        arg1 = varargin{ri+2};
        arg2 = varargin{ri+3};
        tol = varargin{ri+4};
        ri = ri + 5;
        
        r = -1;
        %parsing end and start points
        if numel(arg1) == 2
            angle_start = atan2(arg1(2), arg1(1));
            r = norm( center - arg1 );
        elseif numel(arg1) == 1
            angle_start = arg1;
        else
            error('Invalid start');
        end
        if numel(arg2) == 2
            angle_end = atan2(arg2(2), arg2(1));
            r = norm(center - arg2);
        elseif numel(arg2) == 1
            angle_end = arg2;
        else
            error('Invalid end');
        end
        if (r<0) && ~isnumeric(varargin{ri})
            error('Invalid arc syntax');
        elseif r<0
            r = tol;
            tol = varargin{ri}; ri = ri + 1;
        end
        
        angleDiff = angleDifference(angle_end, angle_start);
        Np = ceil( r*abs(angleDiff) / tol ) + 1;
        angles = angle_start + linspace(0, angleDiff, Np);
        angles = angles(1:(end-1));
        addx( r*[cos(angles); sin(angles); 0*angles] );
        
        if strcmpi(varargin{ri}, 'linename')
            addl(1:(Np-1), varargin{ri+1});
            ri = ri + 2;
        end   
    else
        break;
    end
end

%varargin{ri:end}

if nln > 0
    this.addSurface(x(:,1:np), varargin{ri}, 'Linenames', lnames{1:nln}, varargin{(ri+1):end});
    %lnames{1:nln}
    %disp('dsdas')
    %varargin{(ri+1):end}
else
    this.addSurface(x(:,1:np), varargin{ri:end});
end
end 