function MR = assemble_AGErotationMatrix(varargin)

if numel(varargin) == 2
    %only assembling one time-instant
    rotorAngle = varargin{1};
    N_pp = varargin{2};
    Nterms = 2*N_pp;

    I = zeros(1, 4*N_pp); J = I; E = I;

    for kt = 1:N_pp
        inds = (((kt-1)*4+1):(kt*4));
        I(inds) = (kt-1)*2+[1 1 2 2];
        J(inds) = (kt-1)*2+[1 2 1 2];

        x0 = kt*rotorAngle;
        E(inds) = [cos(x0) -sin(x0) sin(x0) cos(x0)];
    end

    MR = blkdiag(speye(Nterms, Nterms), sparse(I, J, E, Nterms, Nterms));
elseif numel(varargin) == 1
    %assembling a series
    N_pp = varargin{1};
    Nterms = 2*N_pp;
    MR = cell(Nterms, 1);
    MR{1} = blkdiag(speye(Nterms, Nterms), sparse(Nterms, Nterms));
    
    for kp = 2:Nterms
        kf = floor( (kp)/2 );
        inds_row = (((kf-1)*2+1):(kf*2));
        if mod(kp-1,2)
            %sine term
            MR{kp} = blkdiag(sparse(Nterms, Nterms), sparse(inds_row, inds_row([2 1]), [-1 1], Nterms, Nterms));
        else
            %cosine term
            MR{kp} = blkdiag(sparse(Nterms, Nterms), sparse(inds_row, inds_row, [1 1], Nterms, Nterms));
        end
    end
else
    error('Incorrect number of input arguments!');
end
end