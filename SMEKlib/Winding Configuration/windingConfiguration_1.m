function W = windingConfiguration_1(q, p, varargin)
%windingConfiguration_1 Returns a typical 3-phase winding configuration.
% 
% W = windingConfiguration_1(q, p) returns a typical winding configuration W; with
% W(slot) = a ==> phase a traverses slot a to the  sign(a) (positive/negative)
% direction
%   
%   windingConfiguration_1(q, p) returns the winding configuration for a
%   single layer 2p-pole winding with q slots per pole and phase
%   
%   windingConfiguration_1(q, p, a, c) returns the winding configuration for a
%   typical 2p-pole double layer winding with a parallel paths, chorded by c slots
%
% Copyright (c) 2015-2016 Antti Lehikoinen / Aalto University

nargin = length(varargin);

if nargin == 0
    indexMatrix = repmat([1 -3 2 -1 3 -2], 1, p);
    W = kron(indexMatrix, ones(1, q));
elseif nargin == 2    
    chordingPitch = varargin{2};
    a = varargin{1};
    
    if a == 1
        indexMatrix = repmat([1 -3 2 -1 3 -2], 2, 1);
        indexMatrix = repmat(indexMatrix, 1, p);
    elseif a == 2
        indexMatrix = [1 -5 3 -2 6 -4;
            2 -6 4 -1 5 -3];
        indexMatrix = repmat(indexMatrix, 1, p);
    elseif ~mod(a/(2*p), 1)
        % assembling elementary pole-pair indexMatrix
        %{
        num_a = 1; num_b = a + 1; num_c = 2*a+1; %numbering of first path of each phase
        indexMatrix_elem = ...
            [num_a -num_c num_b -num_a-1 num_c+1 -num_b-1;
            num_a+1 -num_c-1 num_b+1 -num_a num_c -num_b];
        indexMatrix = kron(0:(p-1), 2*sign(indexMatrix_elem)) + repmat(indexMatrix_elem, 1, p);
        %}
        indexMatrix = zeros(2, 6*p);
        startSlots = [1 3 5];
        for k_phase = 1:3
            slot_prev = startSlots(k_phase); sign_prev = 1;
            for k_path = 1:a
                indexMatrix(1, slot_prev) = sign_prev*((k_phase-1)*a + k_path);
                slot_prev = mod(slot_prev + 3 - 1, 6*p) + 1;
                indexMatrix(2, slot_prev) = -sign_prev*((k_phase-1)*a + k_path);
                sign_prev = -1*sign_prev;
            end
        end 
    else
        error('The number of parallel paths and/or phases not yet implemented, lol.')
    end

    W = kron(indexMatrix, ones(1, q));

    columnIndices = 1:(size(W,2));
    columnIndices = mod(columnIndices-1+chordingPitch, size(W,2)) + 1;
    W(2,:) = W(2, columnIndices);
else
    error('Incorrect number of input arguments.');
end

end