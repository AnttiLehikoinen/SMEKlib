function DL = doubleLayerFunction_1(slotsPerPoleAndPhase, chordingPitch)

%original
%%{
indexMatrix = [1 -5 3 -2 6 -4;
    2 -6 4 -1 5 -3];
%}

%phases b and c reversed
%{
indexMatrix = [1 5 -3 -2 -6 4;
    2 6 -4 -1 -5 3];
%}

%phase b reversed
%{
indexMatrix = [1 5 3 -2 -6 -4;
    2 6 4 -1 -5 -3];
%}

%phases b and c swapped
%{
indexMatrix = [1 -3 5 -2 4 -6;
    2 -4 6 -1 3 -5];
%}




DL = kron(indexMatrix, ones(1, slotsPerPoleAndPhase));

columnIndices = 1:(size(DL,2));
columnIndices = mod(columnIndices-1+chordingPitch, size(DL,2)) + 1;
DL(2,:) = DL(2, columnIndices);

end