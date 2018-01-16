function pcoeffs = get_NedelecPolynomialCoefficients(elementType, shapeFunction)
%get_NedelecPolynomialCoefficients returns the vector-valued polynomials
%spanning the Nedelec space. NOTE: not the shape functions!
% 
% (c) 2017 Antti Lehikoinen / Aalto University

if elementType == Elements.prism
    %from the book Finite Element Methods for Maxwell's Equations, pp. 200
    if shapeFunction == ShapeFunctions.nedelec
        pcoeffs = cell(1, 9);
        [pcoeffs{:}] = deal(zeros(3,6));
        
        %constants
        pcoeffs{1}(1,1) = 1; %a1
        pcoeffs{5}(2,1) = 1; %a5
        pcoeffs{7}(3,1) = 1; %a7
        
        %x1 i.e. x
        pcoeffs{3}(2,2) = -1; %a3
        pcoeffs{8}(3,2) = 1; %a8
        
        %x2 i.e. y
        pcoeffs{3}(1,3) = 1; %a3        
        pcoeffs{9}(3,3) = 1; %a9
        
        %x3 i.e. z
        pcoeffs{2}(1,4) = 1; %a2
        pcoeffs{6}(2,4) = 1; %a6
        
        %x1*x3, i.e. x*y
        pcoeffs{4}(2,5) = -1; %a4
        
        %x2*x3, i.e. y*z
        pcoeffs{4}(1,6) = 1; %a4   
    else
        error('Element order not yet implemented.');
    end
else
    error('Element shape not yet implemented.');
end

end
       