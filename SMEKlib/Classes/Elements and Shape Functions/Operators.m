classdef Operators < handle
	%{
    enumeration
        %operator types
        div, grad, curl, I
    end
	%}
    methods (Static)		
		%ugly non-enumeration-workaround for enumerations
		function e = div(); e = 20; end;
		function e = grad(); e = 21; end;
		function e = curl(); e = 22; end;
		function e = I(); e = 23; end;
	end
end