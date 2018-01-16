classdef IDfun < handle
    properties
        k
    end
    methods
        function this = IDfun(varargin)
            if numel(varargin)
                this.k = varargin{1};
            else
                this.k = 1;
            end
        end
        function N = eval(~, varargin)
            N = 1;
        end
        function [Nf, order, Nvars] = getData(~, ~)
            Nf = 1; order = 0; Nvars = 1;
        end
        function inds = getIndices(this, ~, ~, ~)
            inds = this.k;
        end
    end
end