function MR = assemble_rotorRotationMatrix(varargin)
% MR = assemble_rotorRotationMatrix(rotorAngle, N_pp)
% assembles the sliding interface rotation matrix for one instant
% 
% MR = assemble_rotorRotationMatrix(N_pp)
% assembles the time-Fourier series of rotation matrices

if numel(varargin) == 2
    %only assembling one time-instant
    rotorAngle = varargin{1};
    N_pp = varargin{2};
    Nterms = 2*N_pp+1;

    I = zeros(1, 1 + 4*N_pp); J = I; E = I;

    I(1) = 1; J(1) = 1; E(1) = 1;
    for kt = 1:N_pp
        inds = 1 + (((kt-1)*4+1):(kt*4));
        I(inds) = 1+(kt-1)*2+[1 1 2 2];
        J(inds) = 1+(kt-1)*2+[1 2 1 2];

        x0 = kt*rotorAngle;
        E(inds) = [cos(x0) -sin(x0) sin(x0) cos(x0)];
    end

    MR = sparse(I, J, E, Nterms, Nterms);
elseif numel(varargin) == 1
    %assembling a series
    N_pp = varargin{1};
    Nterms = 2*N_pp+1;
    MR = cell(Nterms, 1);
    MR{1} = sparse(1, 1, 1, Nterms, Nterms);
    
    for kp = 2:Nterms
        kf = floor( (kp)/2 );
        inds_row = 1 + (((kf-1)*2+1):(kf*2));
        if mod(kp-1,2)
            %sine term
            MR{kp} = sparse(inds_row, inds_row([2 1]), [1 -1], Nterms, Nterms);
        else
            %cosine term
            MR{kp} = sparse(inds_row, inds_row, [1 1], Nterms, Nterms);
        end
    end
else
    error('Incorrect number of input arguments!');
end

end