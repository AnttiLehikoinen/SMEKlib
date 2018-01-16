function this = gw_addSurface(this, points, surfacename, varargin)

%checking if tolerance given
if numel(varargin) && isnumeric(varargin{1})
    [points, ref_inds] = gw_addTolerances(points, varargin{1});
    other_args = {varargin{2:end}};
else
    other_args = varargin;
    ref_inds = [];
end

%adding points to the wrapper object
Np = size(points, 2);
Np_orig = this.N_points;    
this.addPoints(points);

%adding lines
ldef = [1:Np;2:Np 1] + Np_orig;
Nl_orig = this.N_lines;
this.addLines(ldef);

%checking if named lines given
if numel(other_args) && strcmpi(other_args{1}, 'linenames')
    line_inds = (Nl_orig+1):(Nl_orig + Np);
    kend = 1;
    for ri = 2:2:numel(other_args)
        if isnumeric(other_args{ri})
            if isempty(ref_inds)
                lineNumbers = other_args{ri};
            else
                lineNumbers = find( ismember(ref_inds, other_args{ri}) );
            end
            this.physicalLines.add(other_args{ri+1}, line_inds(lineNumbers));
            kend = kend + 2;
        else
            break;
        end
    end
    other_args = {other_args{kend:end}};
end

%adding line loop
lldef = (1:Np) + Nl_orig;
this.addLineloop(lldef);

%adding surface
if numel(other_args) && islogical(other_args{1})
    % surface is a hole
    this.surfaces.add(surfacename, -this.N_lineloops);
else
    this.surfaces.add(surfacename, this.N_lineloops);
    %adding internal surfaces
    sname_ll = ['ll_' num2str(this.N_lineloops)];
    for kh = 1:numel(other_args)
        ntemp = other_args{kh};
        if ischar(ntemp)
            this.holesInSurfaces.add(sname_ll, this.surfaces.get(ntemp));
        else
            this.holesInSurfaces.add(sname_ll, ntemp);
        end
    end
end

this.N_surfaces = this.N_surfaces + 1;

end




