function [hconst, varargout] = lockstepConvolution_2(Xs, h_imp, Nrows, N1, N2, varargin)
%lockstepConvolution_2 convolution with equal time-step lengths.
%
% Vectorized evaluation of several matrix-to-vector convolutions
% simultaneously.
% 
% Each convolution is of 2-input-2-output type; typically from boundary
% vector potentials a and currents i to the Lagrange multiplies h and
% voltages u.
%
% The output vector is of type
% [h1 h2 ... hn; u1 u2 ... un].
%
% It is assumed that each map : (a,i)-->(h,u) is modelled by the same
% matrix-valued discrete impulse response function h_imp; with each row
% containing one time-step value of the vectorized matrix 
% reshape([H;U], [], 1).
%
% The vector Xs contains the inputs; each column corresponding to a single
% time-step with the variables ordered as
% [a1; a2; ... ; an; u1; ... ; un].
%
% Nrows = size(h1,1) + size(u1,1)
% N1 = size(a1,1)
% N2 = size(u1,1)

N_imp = size(h_imp, 1);
N_source = size(Xs,2);
N = N1 + N2;

%determining indices
if nargout > 1
    inds_const_X = max(N_source - N_imp + 1, 1):(N_source-1);
    inds_con st_h = (2 + numel(inds_const_X) - 1):-1:2;
else
    inds_const_X = max(N_source - N_imp + 1, 1):N_source;
    inds_const_h = (1 + numel(inds_const_X) - 1):-1:1;
end

%constant term
N_c = numel(inds_const_X);
Nb = size(Xs,1) / (N1+N2);

%hconst = reshape( transpose(h_imp(inds_const_h,:)), Nh, N_c*N ) * reshape( Xs(:,inds_const_X), N_c*N, 1);
hconst = reshape( transpose(h_imp(inds_const_h,:)), Nrows, N_c*N ) * ...
    reshape(permute(...
    [reshape(transpose(Xs(1:(Nb*N1),inds_const_X)), N_c, N1, Nb) reshape(transpose(Xs((Nb*N1+1):end,inds_const_X)), N_c, N2, Nb)]...
    , [2 1 3]), N_c*N, Nb);

%hconst = reshape(hconst, [], 1);


%including decay term if supplied and not yet fully decayed
if (numel(varargin) > 0) && (size(varargin{1},2) >= N_source)
    %hconst = hconst + varargin{1}(:, N_source);
    hconst = hconst + reshape(varargin{1}(:, N_source), Nrows, Nb);
end

if nargout > 1
    varargout{1} = reshape( transpose(h_imp(1,:)), Nrows, [] );
end

end