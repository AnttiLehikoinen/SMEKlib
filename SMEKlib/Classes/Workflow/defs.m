classdef defs < handle
    %{
	enumeration
        stranded, solid, cage, solidConductor,
        asynchronous, synchronous
        timeharmonic, timestepping
        star, delta
    end
	%}
	methods (Static)		
		%ugly non-enumeration-workaround for enumerations
		function e = stranded();  e = 30; end;
		function e = solid(); e = 31; end;
		function e = cage(); e = 32; end;
		function e = asynchronous(); e = 33; end;
		function e = synchronous(); e = 34; end;
		function e = timeharmonic(); e = 35; end;
		function e = timestepping(); e = 36; end;
		function e = star(); e = 37; end;
		function e = delta(); e = 38; end;
	end
end