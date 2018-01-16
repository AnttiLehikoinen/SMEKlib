function d = matrixDeterminant(M)
%matrixDeterminant3D Matrix determinant (3x3).
% 
% matrixDeterminant3D(M)
%
% Referenced from https://en.wikipedia.org/wiki/Determinant
% 
% (c) 2017 Antti Lehikoinen / Aalto University

if size(M,1) == 4
    d = M(1,:).*M(4,:) - M(2,:).*M(3,:);
elseif size(M,1) == 9
    inds = reshape(1:9, 3, 3);

    d = M(inds(1,1),:).*(M(inds(2,2),:).*M(inds(3,3),:) - M(inds(3,2),:).*M(inds(2,3),:)) - ...
        M(inds(1,2),:).*(M(inds(2,1),:).*M(inds(3,3),:) - M(inds(3,1),:).*M(inds(2,3),:)) + ...
        M(inds(1,3),:).*(M(inds(2,1),:).*M(inds(3,2),:) - M(inds(3,1),:).*M(inds(2,2),:));
end
    

end