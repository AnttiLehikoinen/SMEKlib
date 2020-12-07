figure(10); clf; hold on; box on; axis equal;

%plotting airgap triangulation
AGT = mshc.bandData;
[tag, pag] = mshc.bandData.t_ag( pi/180 * 10 );
triplot(tag', pag(1,:), pag(2,:), 'm');