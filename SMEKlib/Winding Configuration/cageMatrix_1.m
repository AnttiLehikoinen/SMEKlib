function Cr = cageMatrix_1(Qr, varargin)
% Cr = cageMatrix_1 Returns the loop matrix for a rotor cage with Qr
% bars
%   Cr = cageMatrix_1(Qr) ==> end-ring impedance assumed zero
%   Cr = cageMatrix_1(Qr, Zer) == end-ring impedance Z_er between two bars

if isempty(varargin)
    Cr = [ones(1,Qr-1); -eye(Qr-1)];
else
    error('HAHA F U NOT IMPLEMENTED YET')
end

end

