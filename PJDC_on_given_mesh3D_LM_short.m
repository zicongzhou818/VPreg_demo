function [phi1_comp,phi2_comp,phi3_comp,phi_m1, phi_m2, phi_m3]=PJDC_on_given_mesh3D_LM_short(Nx,Ny,Nz,given_x,given_y,given_z,target_x,target_y,target_z)
%% initialization
better = 1;
ei = 0;
ti = 0;
imax = 2.0e2;
ts = 1e-3;
ts_r = 5e-5;
ts_up = 1.1;
ts_dn = 0.9;

r = 1;
r_tol = 1e-8;

%% preset new deformation based on control functions that are corresponsive to the identity map
[X,Y,Z] = ndgrid(1:Nx,1:Ny,1:Nz);
phi_m1 = X;
phi_m2 = Y;
phi_m3 = Z;
phi_m1n = X;
phi_m2n = Y;
phi_m3n = Z;

%% reading given mesh 
phi_o1=given_x;
phi_o2=given_y;
phi_o3=given_z;

%% (ps: the final phi initially is the same as given phi here because the current phi is id)
phi1_comp=given_x;
phi2_comp=given_y;
phi3_comp=given_z;

%% compute ssd
ssd_initial=sum((phi1_comp-target_x).^2+(phi2_comp-target_y).^2+(phi3_comp-target_z).^2,"all");
ssd_old=ssd_initial;

%% main loop
while ti<imax && ts>ts_r && r>r_tol
    ti=ti+1;
    if better % form gradients wrt F when ssd decreases
        ei=ei+1;
        %% form B
        [Gx,Gy,Gz] = computing_G_phi(phi1_comp, phi2_comp, phi3_comp, phi_m1, phi_m2, phi_m3, phi_o1, phi_o2, phi_o3, target_x, target_y, target_z); 
    end
    phi_m1n(2:Nx-1,2:Ny-1,2:Nz-1) = phi_m1(2:Nx-1,2:Ny-1,2:Nz-1) - ts * Gx(2:Nx-1,2:Ny-1,2:Nz-1);
    phi_m2n(2:Nx-1,2:Ny-1,2:Nz-1) = phi_m2(2:Nx-1,2:Ny-1,2:Nz-1) - ts * Gy(2:Nx-1,2:Ny-1,2:Nz-1);
    phi_m3n(2:Nx-1,2:Ny-1,2:Nz-1) = phi_m3(2:Nx-1,2:Ny-1,2:Nz-1) - ts * Gz(2:Nx-1,2:Ny-1,2:Nz-1);
    phi1_compn = interpn(X, Y, Z, phi_o1, phi_m1n, phi_m2n, phi_m3n, 'makima'); 
    phi2_compn = interpn(X, Y, Z, phi_o2, phi_m1n, phi_m2n, phi_m3n, 'makima'); 
    phi3_compn = interpn(X, Y, Z, phi_o3, phi_m1n, phi_m2n, phi_m3n, 'makima'); 
    ssd_new = sum((phi1_compn-target_x).^2 + (phi2_compn-target_y).^2 + (phi3_compn-target_z).^2,"all");
    r = ssd_new/ssd_initial;
    if ssd_new < ssd_old
        ts = ts * ts_up;
        phi_m1 = phi_m1n;
        phi_m2 = phi_m2n;
        phi_m3 = phi_m3n;
        phi1_comp = phi1_compn;
        phi2_comp = phi2_compn;
        phi3_comp = phi3_compn;
        ssd_old = ssd_new;
        better = 1;
    else
        better = 0;
        ts = ts * ts_dn;
    end 
end
display(['Phi_inv found ts:',num2str(ts),' r:',num2str(ssd_new/ssd_initial), ' ei/ti:',num2str(ei), '/',num2str(ti)]);
end
