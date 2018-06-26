function b = matrixTimesVector(M, v, trps, inverse, varargin)

if size(M,1) == 1
    if inverse
        b = v./M;
    else
        b = M.*v;
    end
elseif size(M,1) == 4
    if trps
        ind_11 = 1;
        ind_21 = 3;
        ind_12 = 2;
        ind_22 = 4;
    else
        ind_11 = 1;
        ind_21 = 2;
        ind_12 = 3;
        ind_22 = 4;
    end
    
    %adding constant vector if needed
    if numel(varargin) >= 2
        v = bsxfun(@plus, v, varargin{2});
    end
    if inverse
        %getting determinant
        if ~numel(varargin) || ~any(varargin{1})
            detF = matrixDeterminant(M);
        else
            detF = varargin{1};
        end

        %returns vectorized det(F)*F\V or det(F)*F'\V
        b = [M(ind_22,:).*v(1,:)-M(ind_12,:).*v(2,:);
            -M(ind_21,:).*v(1,:)+M(ind_11,:).*v(2,:)];
        b = bsxfun(@rdivide, b, detF);
    else
        %returns vectorized F*V or F'*B
        b = [M(ind_11,:).*v(1,:)+M(ind_12,:).*v(2,:);
            M(ind_21,:).*v(1,:)+M(ind_22,:).*v(2,:)];
    end    
    
elseif size(M,1) == 9
    %3D mapping

    inds = reshape(1:9, 3, 3);
    if trps
        inds = inds';
    end

    if inverse
        %pain cometh here
        % referenced from https://en.wikipedia.org/wiki/Invertible_matrix
        row1 = crossProduct( M(inds(:,2),:), M(inds(:,3), :) );
        row2 = crossProduct( M(inds(:,3),:), M(inds(:,1), :) );
        row3 = crossProduct( M(inds(:,1),:), M(inds(:,2), :) );
        if ~numel(varargin) || ~any(varargin{1})
            detF = matrixDeterminant(M);
        else
            detF = varargin{1};
        end
        if size(v,2) == 1
            b = [transpose(v)*row1;
                transpose(v)*row2;
                transpose(v)*row3];
        else
            b = [sum(row1.*v,1);
                sum(row2.*v,1);
                sum(row3.*v,1)];
        end
        b = bsxfun(@times, b, 1./detF);
        %b = bsxfun(@rdivide, b, detF);    
    else
        b = bsxfun(@times, M(inds(:,1),:), v(1,:)) + ...
            bsxfun(@times, M(inds(:,2),:), v(2,:)) + ...
            bsxfun(@times, M(inds(:,3),:), v(3,:));
    end
end

end