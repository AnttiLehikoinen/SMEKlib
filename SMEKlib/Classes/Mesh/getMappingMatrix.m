function [F, F0] = getMappingMatrix(this, varargin)

%FIXME make this into an independent function instead
if this.elementType == Elements.triangle || this.elementType == Elements.triangle2
    %getMappingMatrix A mesh-specific mapping matrix.            
    [F, F0] = mappingTerms(this, varargin{:});
elseif this.elementType == Elements.triangle2I
    % [dF, F0] = getMappingMatrix(elements, x0)
    %isoparametric mapping F; returning dF/dxref(x0) and F(x0)
    Nref = Nodal2D(); gradN = Nodal2D(Operators.grad);
    if isempty(varargin{1})
        elements = 1:size(this.t,2);
    else
        elements = varargin{1};
    end
    x0 = varargin{2};
    F = zeros(4, numel(elements)); F0 = zeros(2, numel(elements));
    for kf = 1:size(this.t,1)
        ngrad = gradN.eval_ref(kf, x0, this);
        F = F + [this.p(:, this.t(kf,elements))*ngrad(1);
            this.p(:, this.t(kf,elements))*ngrad(2)];
        F0 = F0 + this.p(:, this.t(kf,elements))*Nref.eval_ref(kf, x0, this);
    end
end

end