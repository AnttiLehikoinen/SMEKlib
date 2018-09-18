Nslices = 1;

angle_skew = 1.5 * 2*pi/dim.Qs;
skew_angles = linspace(-angle_skew/2, angle_skew/2, Nslices);


%updating whatever's changed
Mur =  -speye(Nu_r)-DRr*Lr*(Zr\Lr')/Nslices;
MJ_slice = -1/dims.leff/Nslices*JF_r; %DC-current density matrix
MU_slice = [1/dt*DRr*transpose(JF_r)/Nslices real(Mur)]; %voltage equations

MJs_slice = -JF_s*Ls; %(stator) current source matrix for the slice
MEs_slice = 1/dt*dims.leff*Ls'*JF_s'/Nslices; %(stator) back-emf matrix for the slice

%constant non-ag stiffness matrix and mass matrix
Scell = cell(Nslices+1); [Scell{:}] = deal(sparse(Np+Nu, Np+Nu));
Mcell = cell(Nslices+1); [Mcell{:}] = deal(sparse(Np+Nu, Np+Nu));
Sagcell = cell(Nslices+1, 1); Jcell = cell(Nslices+1,1);
for kslice = 1:Nslices
    Scell{kslice,kslice} = [sparse(Np, Np) MJ_slice;
        sparse(Nu, Np) real(Mur)];
    Scell{kslice,end} = [MJs_slice; sparse(Nu, Ni)];
    Scell{end,kslice} = sparse(Ni, Np+Nu);
    
    Mcell{kslice,kslice} = 1/dt*[M sparse(Np, Nu);
        DRr*transpose(JF_r)/Nslices sparse(Nu,Nu)];
    Mcell{kslice, end} = sparse(Np+Nu, Ni);
    Mcell{end, kslice} = [MEs_slice sparse(Ni,Nu)];
end
Scell{end,end} = Z; 
Mcell{end,end} = sparse(Ni,Ni);
Sagcell{end} = sparse(Ni,Ni);
Jcell{end} = sparse(Ni,Ni);
Stot = cell2mat(Scell); Mtot = cell2mat(Mcell);


PT = assemble_TotalMasterSlaveMatrix(Np + Nu, P_data, []);
Ptot = blkdiag(...
    kron(eye(Nslices), PT), eye(Ni));

Xsamples_slice = zeros((Np+Nu)*Nslices + Ni, Nsamples);
Xsamples_slice(:, 1) = [repmat(Xtot(1:(Np+Nu)), Nslices, 1);
    Xtot((Np+Nu+1):(Np+Nu+Ni))];
Xfree = zeros(size(Ptot,2), 1);

FL = zeros((Np+Nu)*Nslices + Ni, 1);
Ntot = (Np+Nu)*Nslices + Ni;
NfreeInSlice = (size(Ptot,2)-Ni)/ Nslices;
indsI = ((Np+Nu)*Nslices+1):((Np+Nu)*Nslices + Ni);

tic
for kt = 2:2%Nsamples
    disp(['Time step ' num2str(kt) '...']);
    
    disp('Assembling slice matrices');
    for kslice = 1:Nslices
        S_ag = get_MovingBandMatrix(wm*tsamples(kt) + skew_angles(kslice), msh, Np+Nu);    
        Sagcell{kslice} = S_ag;
    end
    FL = Mtot*Xsamples_slice(:,kt-1) + [zeros((Np+Nu)*Nslices,1); Ufun(tsamples(kt))];
    Qconst = Stot + Mtot + blkdiag(Sagcell{:});
    
    %Newton iteration
    for kiter = 1:2%15
        res_t = zeros(Ntot, 1);
        
        %assembling slice matrices
        for kslice = 1:Nslices
            inds_slice = (1:(Np+Nu)) + (kslice-1)*(Np+Nu);
            inds_slice_free = (1:NfreeInSlice) + (kslice-1)*NfreeInSlice;
            [J, res] = assemble_Jacobian(nu_fun, PT*Xfree(inds_slice_free), [], msh);
            Jcell{kslice} = J;
            res_t(inds_slice) = res;
        end
        
        %computing final residual
        restot = Ptot'*(res_t + Qconst*Ptot*Xfree - FL);
        
        %checking convergence
        resNorm = norm(restot) / norm(FL);
        disp([9 'Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']);
        if resNorm < 1e-6
            break;
        end
        %solving
        Jtot = Ptot'*(blkdiag(Jcell{:}) + Qconst)*Ptot;
        
        dX = Jtot\restot;
        %norm(dX)
        
        Xfree = Xfree - dX;
    end
    Xsamples_slice(:, kt) = Ptot*Xfree;
            
end
toc

Isliced = Xsamples_slice(indsI,:);

%plotting currents (and the currents from harmonic approximation
figure(5); clf; hold on; box on;
hs = plot( tsamples*1e3, 2*Xsamples((end-2):end, :)', 'b' );
hm = plot( tsamples*1e3, 2*Isliced', 'r' );
xlabel('Time (ms)')
ylabel('Phase current (A)')
legend([hs(1) hm(1)], 'Single-slice', '11-slice')