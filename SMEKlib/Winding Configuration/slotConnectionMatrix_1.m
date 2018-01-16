function Ck = slotConnectionMatrix_1(N_conductors, N_series, W, k, varargin)
%slotConnectionMatrix_1 loop matrix for a single slot.
% 
% C_k = slotConnectionMatrix_1(N_conductors, N_series, W, k)
%   
% The function returns the connection/loop matrix for the slot k, having N_conductors
% in total, N_series of which are connected in series, based on the winding
% connection matrix W.
%
% With C_k, the currents flowing in the conductors of slot k can be
% obtained as I_slot = C_k * I_free
% 
% Use slotConnectionMatrix_1(N_conductors, N_series, W, k, k_tp)
% to transpose strands from slot to slot in different ways.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%matrix indicating which "phases" go through this slot
indicatorMatrix = zeros(size(W,1), max(max(abs(W))) );
for row = 1:size(W,1)
    indicatorMatrix(row, abs(W(row,k))) = sign(W(row,k));
end

N_inParallel = N_conductors  / N_series / size(indicatorMatrix,1); %number of strands belonging to one current loop

%elementary connection matrix between strands and one current path
elementaryBlock = kron(ones(N_series, 1), eye(N_inParallel));

%"worst-case" configuration
if numel(varargin) == 0
    Ck = kron(indicatorMatrix, elementaryBlock);
    return;
end

transp = varargin{1};

if transp == 0
    Ck = kron(indicatorMatrix, elementaryBlock);
elseif transp == 1
    % return paths twisted
    elementaryBlock_return = kron(ones(N_series, 1), fliplr(eye(N_inParallel)));
    Ck = kron((indicatorMatrix>0), elementaryBlock) + kron(-(indicatorMatrix<0), elementaryBlock_return);
elseif transp == 2
    % block transposition
    elementaryBlock_return = kron(ones(N_series, 1), ...
        kron([0 1;1 0], eye(N_inParallel/2)) );
    Ck = kron((indicatorMatrix>0), elementaryBlock) + kron(-(indicatorMatrix<0), elementaryBlock_return);
    return;
elseif transp == 3
    % alternating transposition
    if mod(k,2)
        elementaryBlock = kron(ones(N_series, 1), fliplr(eye(N_inParallel)));
    end
    Ck = kron(indicatorMatrix, elementaryBlock);
elseif transp == 4
    % alternating block transposition
    if mod(k,2)
        blockSize = varargin{2};
        elementaryBlock = kron(ones(N_series, 1), ...
            kron(fliplr(eye(blockSize)), eye(N_inParallel/blockSize)) );
    end
    Ck = kron(indicatorMatrix, elementaryBlock);
elseif transp == 5
    % alternating reversed block transposition
    if mod(k,2)
        blockSize = varargin{2};
        elementaryBlock = kron(ones(N_series, 1), ...
            kron(eye(blockSize), fliplr(eye(N_inParallel/blockSize))) );
    end
    Ck = kron(indicatorMatrix, elementaryBlock);
else
    error('Transposition type not implemented.')
end

end