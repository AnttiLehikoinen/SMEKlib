function [linehandles, linecoordinates] = drawFluxLines(msh, A, numberOfLines, varargin)
%drawFluxlines Draws flux lines.
%
% drawFluxlines(msh, A, Nl, args) draws a total of Nl flux lines with the
% "plot" function and input arguments args
%
% If the mesh struct msh has a field "rotel", the call syntax can be
% drawFluxlines(msh, A, Nl, rotorAngle, args)
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

N = numberOfLines;

[p_plot, plotArgs] = aux_Plotting_parseInput(msh, varargin{:});

if size(msh.t,1) > 3
    t2 = aux_mesh(msh);
    msh = struct('t', t2, 'p', msh.p);
end

edgeLengths = [p_plot(:,msh.t(2,:))-p_plot(:,msh.t(1,:)) ...
    p_plot(:,msh.t(3,:))-p_plot(:,msh.t(2,:)) ...
    p_plot(:,msh.t(1,:))-p_plot(:,msh.t(3,:))];

delta = sqrt( min(sum(edgeLengths.^2,1)) ) / numberOfLines / 20; %small number for determining unique points

linehandles = cell(N,1);
linecoordinates = cell(N,1);

%potential values to be plotted
relevantNodes = unique(msh.t(:));
potentials = linspace(min(A(relevantNodes)), max(A(relevantNodes)), N+2);
potentials = potentials(2:(N+1));

hold on;

Nprint = min(100, N); printratio = Nprint / N;
fprintf(1,'%s\n\n',repmat('.',1,Nprint)); %for printing a progress "bar"

for kpot = 1:N
    
    if rand < printratio; fprintf(1,'\b.\n'); end
    
    pot = potentials(kpot);
    lineStruct = struct('X', zeros(2,100), 'Nval', 0);
    
    for elem = 1:size(msh.t,2)    
        %check if the equipotential lines goes through this element at all
        elementNodes = msh.t(:,elem);
        if (max(A(elementNodes)) < pot) || (min(A(elementNodes)) > pot)
            %nope --> continue to the next element in list
            continue;
        end
        
        %drawing the equipotential line in the reference element:
            %denoting nodal values
            a1 = A(elementNodes(1)); a2 = A(elementNodes(2)); a3 = A(elementNodes(3));

            %finding the intersects of the eq. line and the ksi-axis, eta-axis
            %and the 1-ksi line
            ksiis = [(pot-a1)/(a2-a1) 0 (pot-a3)/(a2-a3)];
            etas = [0 (pot-a1)/(a3-a1) 1-(pot-a3)/(a2-a3)];

            %finding the intersection points that are actually inside the
            %triangle
            I = find( (ksiis>=0).*(ksiis<=1).*(etas<=1).*(etas>=0) );
            ksiis = ksiis(I);
            etas = etas(I);
            localCoordinates = [ksiis;etas];
            
            %orienting the eqpot line
            gradA = a1*[-1;-1] + a2*[1;0] + a3*[0;1];
            PV = localCoordinates(:,2)-localCoordinates(:,1);
            if ( PV(1)*gradA(2) - PV(2)*gradA(1) ) < 0
                localCoordinates = localCoordinates(:,[2 1]);
            end
            
            
        %calculating mapping from reference element to global element
        AM = [p_plot(:,elementNodes(2))-p_plot(:,elementNodes(1)) p_plot(:,elementNodes(3))-p_plot(:,elementNodes(1))];
        b = p_plot(:,elementNodes(1)*[1 1]);
        
        
        %A = [msh.x(elementNodes(2))-msh.x(elementNodes(1)) msh.x(elementNodes(3))-msh.x(elementNodes(1));
        %    msh.y(elementNodes(2))-msh.y(elementNodes(1)) msh.y(elementNodes(3))-msh.y(elementNodes(1))];
        %b = [msh.x(elementNodes(1)); msh.y(elementNodes(1))];
        
        %calculating and plotting global equipotential lines
        globalCoordinates = AM*localCoordinates + b;
        
        %plot(globalCoordinates(1,:), globalCoordinates(2,:), style);
        lineStruct = appendValue(globalCoordinates, lineStruct);
    end

    X = lineStruct.X(:, 1:lineStruct.Nval);
    eqLines = lineReorderer(X, delta, numberOfLines);
    linehandles{kpot} = cell(size(eqLines,1), 1); 
    for kl = 1:size(eqLines,1)
        linehandles{kpot}{kl} = plot( eqLines{kl}(1,:), eqLines{kl}(2,:), plotArgs{:});
    end
    
    linecoordinates{kpot} = eqLines;
    %linehandles(kpot) = plot( X(1,:), X(2,:), style);
end
   
%linecoordinates = 0;
%linehandles = 0;

end

function S = appendValue(X, S)

if size(X, 2) ~= 2
    error('Hurrrrr')
end

if (S.Nval+2) <= size(S.X, 2)
    %
else
    S.X = [S.X zeros(2, 3*size(S.X,2))];
end
S.X(:, S.Nval+[1 2]) = X;
S.Nval = S.Nval + 2;

end

function eqLines = lineReorderer(X, delta, numberOfLines)
%getting nodes that are unique within precision 'delta'
[~,IA,IC] = unique(floor(X/delta)', 'rows', 'stable');

X2 = X(:,IA);

lineNodes = [1:2:size(X,2); 2:2:size(X,2)];
lineNodes = [IC(lineNodes(1,:))'; IC(lineNodes(2,:))']; %lines expressed with unique nodes

visitedNodes = zeros(1, numel(IA));

%node connection matrix
A = sparse([lineNodes(1,:) lineNodes(2,:)], [lineNodes(2,:) lineNodes(1,:)], ones(1,2*size(lineNodes,2)), numel(IA), numel(IA));

%going through un-closed lines
orphanNodes = find( sum(A, 2) < 2 );
paths = cell(numberOfLines,1); mainPathIndex = 1;
while 1
    if ~any(orphanNodes)
        break;
    end
    currentNode = orphanNodes(1);
    thisPath = zeros(1, numel(IA)); pathIndex = 1;
    
    while 1
        visitedNodes(currentNode) = 1; 
        thisPath(pathIndex) = currentNode;
        nextNodes = setdiff([find(A(currentNode,:)) find(A(:,currentNode))'], currentNode);
        if ~any(nextNodes) || ~any(~visitedNodes(nextNodes))
            break;
        end
        for k = nextNodes
            if ~visitedNodes(k)
                currentNode = k;
                pathIndex = pathIndex + 1;
                break;
            end
        end
    end
    %disp(thisPath(1:pathIndex))
    paths{mainPathIndex} = thisPath(1:pathIndex); mainPathIndex = mainPathIndex + 1;
    
    orphanNodes = setdiff(orphanNodes, thisPath(1:pathIndex));
end

%going through closed lines
while 1
    currentNode = find( ~visitedNodes, 1);
    if ~any(currentNode)
        break;
    end
    thisPath = zeros(1, numel(IA)); pathIndex = 1;
    
    while 1
        visitedNodes(currentNode) = 1; 
        thisPath(pathIndex) = currentNode;
        nextNodes = setdiff([find(A(currentNode,:)) find(A(:,currentNode))'], currentNode);
        if ~any(nextNodes) || ~any(~visitedNodes(nextNodes))
            break;
        end
        for k = nextNodes
            if ~visitedNodes(k)
                currentNode = k;
                pathIndex = pathIndex + 1;
                break;
            end
        end
    end
    thisPath(pathIndex) = thisPath(1); %pathIndex = pathIndex + 1; %closing path
    
    paths{mainPathIndex} = thisPath(1:pathIndex); 
    mainPathIndex = mainPathIndex + 1;
end


eqLines = cell(mainPathIndex-1, 1);
for kp = 1:(mainPathIndex-1)
    eqLines{kp} = X2(:, paths{kp});
end




end