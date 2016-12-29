function Cs = statorConnectionMatrix(W, N_conductorsPerSlot, N_series, varargin)
%statorConnectionMatrix returns the stator loop matrix.
% 
% Basic call syntax
% Cs = statorConnectionMatrix(W, N_conductorsPerSlot, N_series)
% where W is the winding configuration matrix. Use
% Cs = statorConnectionMatrix(W, N_conductorsPerSlot, N_series, transp)
% to specify transposition of conductors between slots. See the function slotConnectionMatrix_1
% for available transposition types.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

Qs = size(W, 2);
N_phases = max(abs(W(:)));
N_layers = size(W,1);
N_inParallel = N_conductorsPerSlot  / N_series / N_layers;

Cs = zeros(Qs*N_conductorsPerSlot, N_inParallel*N_phases);

for kslot = 1:Qs
    Cs( ((kslot-1)*N_conductorsPerSlot+1):(kslot*N_conductorsPerSlot), : ) = ...
        slotConnectionMatrix_1(N_conductorsPerSlot, N_series, W, kslot, varargin{:}) ;
end

end

