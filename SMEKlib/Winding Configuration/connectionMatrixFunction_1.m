function Ck = connectionMatrixFunction_1(N_strands, N_series, Double_Layer, slot)
%makes the connection matrix for slot nr. "slot"
%a double layer with 2 parallel paths assumed

%matrix indicating which "phases" go through this slot
indicatorMatrix = zeros(size(Double_Layer,1), max(max(abs(Double_Layer))) );
for row = 1:size(Double_Layer,1)
    indicatorMatrix(row, abs(Double_Layer(row,slot))) = sign(Double_Layer(row,slot));
end

%indicatorMatrix = zeros(2, 6);
%indicatorMatrix(1, abs(Double_Layer(slot,1))) = sign(Double_Layer(slot,1));
%indicatorMatrix(2, abs(Double_Layer(slot,2))) = sign(Double_Layer(slot,2));


strandsPerLoop = N_strands  / N_series / size(indicatorMatrix,1); %number of strands belonging to one current loop

%elementary connection matrix between strands and one current path
elementaryBlock = kron(ones(N_series, 1), eye(strandsPerLoop));

%"worst-case" configuration
%Ck = kron(indicatorMatrix, elementaryBlock);

%"best-case" configuration (return path twisted)
%%{
elementaryBlock_return = kron(ones(N_series, 1), fliplr(eye(strandsPerLoop)));
Ck = kron((indicatorMatrix>0), elementaryBlock) + kron(-(indicatorMatrix<0), elementaryBlock_return);
%}

%only upper layer return path twisted
%{
elementaryBlock_return = kron(ones(N_series, 1), fliplr(eye(strandsPerLoop)));
Ck = kron( [(indicatorMatrix(1,:)>0); indicatorMatrix(2,:)], elementaryBlock) + ...
    kron( [-(indicatorMatrix(1,:)<0); 0*indicatorMatrix(2,:)], elementaryBlock_return);
%}

%twisting all upper layer paths
%{
elementaryBlock_return = kron(ones(N_series, 1), fliplr(eye(strandsPerLoop)));
Ck = kron( [indicatorMatrix(1,:);0*indicatorMatrix(2,:)], elementaryBlock_return) + ...
    kron( [0*indicatorMatrix(1,:);indicatorMatrix(2,:)], elementaryBlock);
%}

%twisting all lower layer paths
%{
elementaryBlock_return = kron(ones(N_series, 1), fliplr(eye(strandsPerLoop)));
Ck = kron( [indicatorMatrix(1,:);0*indicatorMatrix(2,:)], elementaryBlock) + ...
    kron( [0*indicatorMatrix(1,:);indicatorMatrix(2,:)], elementaryBlock_return);
%}

end