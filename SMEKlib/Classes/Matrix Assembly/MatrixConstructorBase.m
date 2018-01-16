classdef MatrixConstructorBase < handle
    properties
        I, J, E, Nvals, Ncoord
        Nrows, Ncols
    end
    methods
        function this = MatrixConstructorBase(varargin)
            if numel(varargin)
                INIT_SIZE = varargin{1};
            else
                INIT_SIZE = 1000;
            end
            this.I = zeros(1, INIT_SIZE);
            this.J = zeros(1, INIT_SIZE);
            this.E = zeros(1, INIT_SIZE);
            this.Nvals = 0;
            this.Ncoord = 0;
            this.Nrows = 0;
            this.Ncols = 0;
        end
        
        function this = addValues(this, E, varargin)
            if numel(varargin)
                pos = varargin{1};
            else
                pos = this.Nvals;
            end
            if (pos + numel(E)) > numel(this.E)
                n_add = max(numel(E), 2*numel(this.E));
                this.E = [this.E zeros(1,n_add)];
            end
            if numel(varargin)
                this.E( (pos+1):(pos+numel(E)) ) = this.E( (pos+1):(pos+numel(E)) ) + toRow(E);
            else
                this.E( (pos+1):(pos+numel(E)) ) = toRow(E);
            end
            this.Nvals = max(this.Nvals, pos+numel(E));
        end
        
        function this = addCoordinates(this, I, J)
            if numel(I) == 1
                I = I*ones(1, numel(J));
            end
            if numel(J) == 1
                J = J*ones(1, numel(I));
            end
            if (this.Ncoord + numel(I)) > numel(this.I)
                n_add = max( numel(I), 2*numel(this.I) );
                this.I = [this.I zeros(1,n_add)];
                this.J = [this.J zeros(1,n_add)];
            end
            this.I( (this.Ncoord+1):(this.Ncoord+numel(I)) ) = toRow(I);
            this.J( (this.Ncoord+1):(this.Ncoord+numel(J)) ) = toRow(J);
            this.Ncoord = this.Ncoord + numel(I);
        end
        
        function S = finalize(this, varargin)
            if this.Nvals ~= this.Ncoord
                error('The numbers of coordinate-pairs and values differ.');
            end
            if ~numel(varargin)
                S = sparse(this.I(1:this.Ncoord), this.J(1:this.Ncoord), ...
                    this.E(1:this.Nvals), this.Nrows, this.Ncols);
                return
            end
            S  = sparse(this.I(1:this.Ncoord), this.J(1:this.Ncoord), ...
                this.E(1:this.Nvals), varargin{:});
        end
        
        function this = resetValues(this)
            this.Nvals = 0;
        end
        function this = resetCoordinates(this)
            this.Ncoord = 0;
        end
    end
end