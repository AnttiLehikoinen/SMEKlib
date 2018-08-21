function [H, dHdB] = matfun_1D(B, nu_fun)

[nu, dnu] = nu_fun([B;0*B]);

H = nu.*B;
dHdB = nu+2*dnu.*B.*B;

end