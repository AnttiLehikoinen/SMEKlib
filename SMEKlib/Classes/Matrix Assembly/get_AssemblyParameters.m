function [x_quad, w_quad, Nrows, Ncols, Nf1, Nf2] = ...
                get_AssemblyParameters(fun_test, fun_shape, msh, varargin)
            
[Nf1, order_test, Nrows] = fun_test.getData(msh);
[Nf2, order_shape, Ncols] = fun_shape.getData(msh);

if numel(varargin)
    order_nonlin = varargin{1}; %extra integration order due to eg. nonlinearity
    if order_shape == 0
        order_nonlin = 0; %nonlinearity has no effect if shape functions are constants
    end
    if numel(varargin)>1
        if ~isempty(varargin{2})
            order_nonlin = order_nonlin + varargin{2};
        end
    end
else
    order_nonlin = 0;
end

if Elements.isTriangle( msh.elementType )
    order_curvature = curvatureOrder( msh.elementType );    
    [x_quad, w_quad] = get_2DtriangleIntegrationPoints(order_test + order_shape + order_nonlin + order_curvature);
    %[x_quad, w_quad] = get_2DtriangleIntegrationPoints(3);
elseif Elements.isTet( msh.elementType )
    order_curvature = 0;
    [x_quad, w_quad] = get_3DtetIntegrationPoints(order_test + order_shape + order_nonlin + order_curvature);
    %[x_quad, w_quad] = get_3DtetIntegrationPoints(2);
    %[x_quad, w_quad] = gaussint_tet(order_test + order_shape + order_nonlin + order_curvature);
    %w_quad = transpose(w_quad);
    %[x_quad, w_quad] = get_3DtetIntegrationPoints(order_test + order_shape + order_nonlin + order_curvature);
elseif Elements.isPrism( msh.elementType )
    order_curvature = 0;
    [x_quad, w_quad] = get_PrismIntegrationPoints(order_test + order_shape + order_nonlin + order_curvature);
else
    error('Not yet implemented');
end

end

function n = curvatureOrder(type)
    % extra integration order due to isoparametric element curvature
    if ~Elements.isIsoparametric(type)
        n = 0; return;
    end
    if type == Elements.triangle2I
        n = 1;
        return;
    end
end