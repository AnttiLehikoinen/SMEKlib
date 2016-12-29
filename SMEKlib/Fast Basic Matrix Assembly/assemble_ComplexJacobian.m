function [Jrr, Jri, Jir, Jii, rr, ri] = assemble_ComplexJacobian(nu_fun, X, ~, msh)
%assemble_ComplexJacobian assembles the Jacobian for a complex problem.
% 
% [Jrr, Jri, Jir, Jii, rr, ri] = assemble_ComplexJacobian(nu_fun, X, elements, msh)
% returns the Jacobian
% [Jrr Jri; Jir Jii]
% for solving the complex problem rr + 1i*ri = 0 with Newton's method.
%
% NOTE: "elements" does not have any effect yet.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

Nvars = size(X, 1) / 2;

%N = size(msh.p, 2); %Jacobian size
N = Nvars; 


%computing B and nu
[~, Breal] = calculate_B(X(1:Nvars), msh);
[~, Bimag] = calculate_B(X((Nvars+1):end), msh);

[nu, dnu] = nu_fun( sum(Breal.^2 + Bimag.^2,1) );

%linearity enforcement
%Breal = 0*Breal; Bimag = 0*Bimag;
%[nu, dnu] = nu_fun( 0*sum(Breal.^2 + Bimag.^2,1) ); dnu = 0*dnu;

Hreal = bsxfun(@times, Breal, nu);
Himag = bsxfun(@times, Bimag, nu);

% assembling:

%drr / dXr
dHdB = [nu+2*dnu.*Breal(1,:).*Breal(1,:);
    2*dnu.*Breal(1,:).*Breal(2,:);
    2*dnu.*Breal(2,:).*Breal(1,:);
    nu+2*dnu.*Breal(2,:).*Breal(2,:)];

[E, I, J, E_res, I_res] = assemble_Jacobian_data(@(B)(deal(Hreal, dHdB)), X(1:Nvars), [], msh, true);

Jrr = sparse(I, J, E, N, N);
rr = sparse(I_res, ones(size(I_res)), E_res, Nvars, 1);

%drr / dXi
dHdB = 2*[dnu.*Breal(1,:).*Bimag(1,:);
    dnu.*Breal(2,:).*Bimag(1,:);
    dnu.*Breal(1,:).*Bimag(2,:);
    dnu.*Breal(2,:).*Bimag(2,:)];
E = assemble_Jacobian_data(@(B)(deal([], dHdB)), X((Nvars+1):(2*Nvars)), [], msh, false);
Jri = sparse(I, J, E, N, N);

%dri / dXr
dHdB = [2*dnu.*Bimag(1,:).*Breal(1,:);
    2*dnu.*Bimag(2,:).*Breal(1,:);
    2*dnu.*Bimag(1,:).*Breal(2,:);
    2*dnu.*Bimag(2,:).*Breal(2,:)];
E = assemble_Jacobian_data(@(B)(deal([], dHdB)), X(1:Nvars), [], msh, false);
Jir = sparse(I, J, E, N, N);

%dri / dXi
dHdB = [nu+2*dnu.*Bimag(1,:).*Bimag(1,:);
    2*dnu.*Bimag(1,:).*Bimag(2,:);
    2*dnu.*Bimag(2,:).*Bimag(1,:);
    nu+2*dnu.*Bimag(2,:).*Bimag(2,:)];
[E, E_res] = assemble_Jacobian_data(@(B)(deal(Himag, dHdB)), X((Nvars+1):(2*Nvars)), [], msh, true);
Jii = sparse(I, J, E, N, N);
ri = sparse(I_res, ones(size(I_res)), E_res, Nvars, 1);

end