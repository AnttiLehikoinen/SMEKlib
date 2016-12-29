function [PT, varargout] = assemble_TotalMasterSlaveMatrix(Ntot, P, varargin)
%assemble_TotalMasterSlaveMatrix elimination matrix for master-slave
%couplings.
% 
% [PT, otherNodes_new] = assemble_TotalMasterSlaveMatrix(Ntot, P, otherNodes)
% generates a Ntot x N mapping matrix, with N <= Ntot
% 
% Each of the matrices in the cell array P{1,:} defines a part of the
% mapping.
% 
% Alternatively P{1,c} can be a 3xN sparse triplet array.
%
% Optionally, P{2,:} and P{3,:} can contain index arrays for
% permuting/shifting the rows and columns of the P-matrices, respectively.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%determining free nodes
N_Ps = size(P, 2); %number of interpolation blocks given

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding entries of interpolation matrices

Pinds = cell(3, N_Ps);
if size(P,1) > 1
    checkIndices = true;
else
    checkIndices = false;
end

for kp = 1:N_Ps
    if issparse(P{1, kp})
        [I, J, E] = find(P{1, kp});
    elseif size(P{1,kp},1) == 3
        [I, J, E] = deal( P{1, kp}(1,:)', P{1, kp}(2,:)', transpose(P{1, kp}(3,:)) );
    else
        error('Incorrect input in P.');
    end
    
    %checking if P_kp needs to be row- or column-scaled
    if checkIndices && ~isempty(P{2,kp})
        Pinds{1,kp} = P{2,kp}(I');
    else
        Pinds{1,kp} = I';
    end
    
    if checkIndices && ~isempty(P{3,kp})
        Pinds{2,kp} = P{3,kp}(J');
    else
        Pinds{2,kp} = J';
    end
    Pinds{3,kp} = E';
end

slaveNodes = unique( horzcat( Pinds{1,:} ) );

[freeNodes, IA] = setdiff(1:Ntot, slaveNodes);

Np_free = numel(freeNodes);
indsOfFree(IA) = 1:Np_free;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%constructing total matrix
PT = [];

%free nodes
PT = sparseAdd(freeNodes, 1:Np_free, ones(1,Np_free), PT);

%slave nodes
for kp = 1:N_Ps
    %checking for possibly zero input
    inds = find( Pinds{3,kp} ~= 0 );
    
    PT = sparseAdd( Pinds{1,kp}(inds), indsOfFree( Pinds{2,kp}(inds) ), Pinds{3,kp}(inds), PT );
    
    %any( ~indsOfFree( Pinds{2,kp}(inds) ) )
end
PT = sparseFinalize(PT, Ntot, Np_free);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%updating node indices
if numel(varargin) == 0;
    return;
end

varargout = cell( size(varargin) );
for k = 1:numel(varargin)
    varargout{k} = indsOfFree( varargin{k} );
end

end