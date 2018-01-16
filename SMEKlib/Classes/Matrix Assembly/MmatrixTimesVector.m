function b = MmatrixTimesVector(M, v, trps, inverse, varargin)

if size(M,1) == 4
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
    
elseif size(M,1) == 6 && size(v,1) == 3
    %purely for magnetomechanical problems; no transposing or inverses
    inds = reshape(1:6, 2, 3);
    b = bsxfun(@times, M(inds(:,1),:), v(1,:)) + ...
                bsxfun(@times, M(inds(:,2),:), v(2,:)) + ...
                bsxfun(@times, M(inds(:,3),:), v(3,:));
elseif size(M,1) == 6 && size(v,1) == 2
    inds = reshape(1:6, 3, 2);
    b = bsxfun(@times, M(inds(:,1),:), v(1,:)) + ...
                bsxfun(@times, M(inds(:,2),:), v(2,:));
elseif size(M,1) == 9
    %purely for magnetomechanical problems; no transposing or inverses
    inds = reshape(1:9, 3, 3);
    b = bsxfun(@times, M(inds(:,1),:), v(1,:)) + ...
                bsxfun(@times, M(inds(:,2),:), v(2,:)) + ...
                bsxfun(@times, M(inds(:,3),:), v(3,:));
end

end