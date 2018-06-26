function L = EWsegmentIncidenceMatrix(W, Nturns)
%EWsegmentIncidenceMatrix Indicence matrix from phase to end-winding
%segment.
%
% L = EWsegmentIncidenceMatrix(Q, W, Nturns), where
%   W = winding configuration matrix
%   Nturns = number of turns per slot and phase
%
% Function returns the incidence matrix L with the entries
% L(r,m) = n, where
%   r = number of end-winding segment = segment between slots r and r+1
%   m = number of phase
%   n = number of turns per ew-segment and phase.
%
% NOTE: only works for two-layer windings, for now.
%
% (c) 2018 Antti Lehikoinen / Aalto University

Q = size(W,2); %number of slots
m = max(W(:)); %number of phases

L = zeros(Q, m);
Wtemp = W;
%going through slots
for kslot = 1:Q
    
    %going through "layers" = layer 2 = bottom layer by assumption
    for klayer = 2
        ind_m = W(klayer, kslot); %phase that goes here
        Wtemp(klayer, kslot) = 0; %this is considered
        
        %searching for return path
        for kslot2 = 1:Q
            %if/when we flip from slot Q to 1 etc
            kslot_actual = mod( kslot2 + kslot -2, Q)+1;
            
            %return path here in layer 1 = top layer
            if Wtemp(1,kslot_actual)==-ind_m
                Wtemp(1, kslot_actual) = 0;
                break;
            end
            L(kslot_actual, abs(ind_m)) = L(kslot_actual, abs(ind_m)) + sign(ind_m);
        end
    end
end
L = Nturns*L;

end