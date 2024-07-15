function [phi, I_m_phi, phi_inv, I_f_phi_inv] = ImgReg_SymVPcontrol(I_m,I_f,J_thr,regrid_bd,name_is)
I_m_orig = I_m;
I_f_orig = I_f;
I_m_max=max(I_m,[],"all");
I_f_max=max(I_f,[],"all");
I_m=I_m/I_m_max;
I_m=I_m*I_f_max;
[Nx,Ny,Nz] = size(I_m);
[Mx,My,Mz] = size(I_f);
sizes = [Nx,Ny,Nz,Mx,My,Mz];
Nmax = max(sizes,[],"all");
Nmin = min(sizes,[],"all");
N_avg = (Nmin+Nmax)/4;
if rem(N_avg,2)==0
    N_avg=N_avg+1;
elseif rem(N_avg,2)~=0
    N_avg=ceil(N_avg);
end
alternating = regrid_bd - 2;
%% Image registration 
stage = 0;
stagnate = 10;
ts_initF = 1e-6;
ts_init = 1e-5;
ts = ts_initF;% tstep: initial tstep 
ts_rtol = 1e-10;% tstep_ratio: stop when tstep less than tstep_ratio
mod=['F', '-mod'];
lr_ctlF = 0.1;
lr_ctl = 0.05;
ts_up = 1 + lr_ctlF;% tstep_up: magnification factor when ssd decreases
ts_dn = 1 - lr_ctlF;% tstep_down: reduction factor when doesn't ssd decrease
imax = 1.00e2;% ite_max: preset maximum number of iterations
rtol = 6e-1;% ratiotol: preset ratio to stop 
r = 1;
% running time by tic and toc
% h = 2;
%% Step-1
grid_size = [Mx-2,My-2,Mz-2];% remove zero boundaries
F = zeros([grid_size, 3]);
f = ones(grid_size);
g = zeros([grid_size, 3]);

% Compute initial ssd on uniform grid
ssd_old = sum((I_m-I_f).^2,"all");
ssd_initial = ssd_old;
ssd_initial_temp = ssd_old;
% ssd_improve = 1;
mi = mutualInfo(I_m, I_f);
% iteration flag
Better = 1;% Better: indicate whether the ssd decreases(Better=1) or not 
true_better = 0;
ti = 0;% ii: total iteration steps
ei = 0;% iii: iteration steps that ssd decreases

% set up initial values for Phi(x,y)=(x+u1(x,y),y+u2(x,y))
[X,Y,Z] = ndgrid(1:Mx,1:My,1:Mz);% (i,j)--->(pos_x(i,j),pos_y(i,j)) on (1 to m,1 to n)
phi_x = X;
phi_y = Y;
phi_z = Z;

phi_x_interm = X;
phi_y_interm = Y;
phi_z_interm = Z;

I_m_local = I_m;
JacD_below = 1;
lll = 0;
tic;
while (lll <= regrid_bd) 
    while  (r > rtol) && (ts > ts_rtol) && (stage < stagnate) && (ei < imax) && (JacD_below >= J_thr) 
        ti=ti+1; % record the number of iterations
        %% Step-2
        if Better  % --(PS: Better holds for the first step in "while" loop)
            %% 3
            ei=ei+1;
            I_m_temp = interpn(X,Y,Z,I_m_local,phi_x, phi_y,phi_z,'makima'); 
            Img_diff = I_m_temp-I_f; 
            [I_m_y, I_m_x, I_m_z] = gradient(I_m_temp); 
            lap_B(:,:,:,1) = Img_diff(2:Mx-1,2:My-1,2:Mz-1).*I_m_x(2:Mx-1,2:My-1,2:Mz-1);
            lap_B(:,:,:,2) = Img_diff(2:Mx-1,2:My-1,2:Mz-1).*I_m_y(2:Mx-1,2:My-1,2:Mz-1);
            lap_B(:,:,:,3) = Img_diff(2:Mx-1,2:My-1,2:Mz-1).*I_m_z(2:Mx-1,2:My-1,2:Mz-1);
            parfor ith=1:3
                    RHS = imresize3(lap_B(:,:,:,ith),[N_avg-2 N_avg-2 N_avg-2]);
                    B(:,:,:,ith) = pois3fft(RHS);                
            end
            if lll<alternating 
                parfor ith=1:3
                    del_F(:,:,:,ith) = imresize3(B(:,:,:,ith),[Mx-2 My-2 Mz-2]);
                end
            else
                parfor ith=1:3
                    del_F(:,:,:,ith) = imresize3(B(:,:,:,ith),[Mx-2 My-2 Mz-2]);
                end
                [B1_x2,  B1_x1,  B1_x3] = gradient(del_F(:,:,:,1)); 
                [B2_x2,  B2_x1,  B2_x3] = gradient(del_F(:,:,:,2)); 
                [B3_x2,  B3_x1,  B3_x3] = gradient(del_F(:,:,:,3)); 
                del_f = - f.*(B1_x1 + B2_x2 + B3_x3);
                del_g(:,:,:,1) = - B3_x2 + B2_x3;
                del_g(:,:,:,2) =   B3_x1 - B1_x3;
                del_g(:,:,:,3) = - B2_x1 + B1_x2;
            end
        end
        if lll<alternating
            F_n = F - del_F * ts;
            parfor ith=1:3
                lap_U(:,:,:,ith) = imresize3(F_n(:,:,:,ith),[N_avg-2 N_avg-2 N_avg-2]);
            end
        else
            f_n = f - 80.0 * del_f * ts;
            g_n = g - 80.0 * del_g * ts;
            %% Form Lap from 
            [f_x2, f_x1, f_x3] = gradient(f_n); 
            [g1_x2, ~, g1_x3] = gradient(g_n(:,:,:,1)); 
            [~, g2_x1, g2_x3] = gradient(g_n(:,:,:,2)); 
            [g3_x2, g3_x1, ~] = gradient(g_n(:,:,:,3)); 
            fg1 = f_x1 - g3_x2 + g2_x3;
            fg2 = f_x2 + g3_x1 - g1_x3;
            fg3 = f_x3 - g2_x1 + g1_x2;
            fg(:,:,:,1) = imgaussfilt3(fg1, 2);
            fg(:,:,:,2) = imgaussfilt3(fg2, 2); 
            fg(:,:,:,3) = imgaussfilt3(fg3, 2); 
            %% Step-4 Solve for U, where lap(U)=F
            parfor ith=1:3
                lap_U(:,:,:,ith) = imresize3(fg(:,:,:,ith),[N_avg-2 N_avg-2 N_avg-2]);
            end
        end
        parfor ith=1:3
            U_is = pois3fft(lap_U(:,:,:,ith));
            U(:,:,:,ith) = matrixpad3D(imresize3(U_is,[Mx-2 My-2 Mz-2]),0);
        end
        %% Step-6  % phi(x,y)=X+U
        phi_x = X + U(:,:,:,1);
        phi_y = Y + U(:,:,:,2);
        phi_z = Z + U(:,:,:,3);

        [JacD2, ~] = compute_JD_and_Curl3D(phi_x, phi_y, phi_z, 1);
        JacD_below = min(JacD2, [], "all");
        I_m_phi = interpn(X, Y, Z, I_m_local, phi_x, phi_y, phi_z,'makima'); 
        mi_fwd = mutualInfo(I_m_phi, I_f);
        mi_inc = (mi_fwd - mi)/mi;
        ssd = sum((I_m_phi-I_f).^2,"all");
        r = ssd/ssd_initial_temp;
        if ssd<ssd_old
            stage = 0;
            I_m_better = I_m_phi;
            ts = ts * ts_up;
            ssd_old = ssd;
            Better = 1;
            if lll<alternating
                F = F_n;
            else
                f = f_n;
                g = g_n;
                mod = ['f&g', '-mod'];
            end
            phi_x_better = phi_x;
            phi_y_better = phi_y;
            phi_z_better = phi_z;
            true_better = 1;
        else
            Better = 0;
            ts = ts*ts_dn;
            stage = stage + 1;
        end 
        % display([mod,' comp#',num2str(lll),' ei/ti:',num2str(ei),'/',num2str(ti),' ts:',num2str(ts),' r/rr:',num2str(ssd/ssd_initial),'/',num2str(ssd/ssd_initial_temp),' J_min:',num2str(JacD_below),' rMI:%', num2str(mi_inc*100)]);
    end
    display([mod,' comp#:',num2str(lll),' ei/ti:',num2str(ei),'/',num2str(ti),' ts:',num2str(ts),' r/rr:',num2str(ssd/ssd_initial),'/',num2str(ssd/ssd_initial_temp),' J_min:',num2str(JacD_below),' rMI:%', num2str(mi_inc*100)]);
    JacD_below = 1;
    if true_better > 0
        I_m_local = I_m_better;
        phi_x_interm = interpn(X, Y, Z, phi_x_interm, phi_x_better, phi_y_better, phi_z_better, 'makima'); 
        phi_y_interm = interpn(X, Y, Z, phi_y_interm, phi_x_better, phi_y_better, phi_z_better, 'makima');
        phi_z_interm = interpn(X, Y, Z, phi_z_interm, phi_x_better, phi_y_better, phi_z_better, 'makima');
        I_m_phi_true = interpn(X, Y, Z, I_m, phi_x_interm, phi_y_interm, phi_z_interm,'makima'); 
        ssd_old = sum((I_m_phi_true-I_f).^2, "all");
        r_fwd = ssd_old/ssd_initial;
        ssd_initial_temp = ssd_old;
        JD_whole = compute_JD3D(phi_x_interm, phi_y_interm, phi_z_interm);
        JD_min_whole = min(JD_whole,[],"all");   
        true_better = 0;
        display([' r_fwd:', num2str(r_fwd),' J_min_comp:', num2str(JD_min_whole)]);
    end
    phi_x = X;
    phi_y = Y;
    phi_z = Z;
    lll = lll + 1;
    stage = 0;
    Better = 1;
    r=1;
    ei=0;
    ti=0;   
    if lll<alternating
        F = zeros([grid_size, 3]);
        ts = ts_initF;
        ts_up = 1 + lr_ctlF;
        ts_dn = 1 - lr_ctlF;
    else
        f = ones(grid_size);
        g = zeros([grid_size, 3]);
        ts = ts_init;
        ts_up = 1 + lr_ctl;
        ts_dn = 1 - lr_ctl;
        rtol = 5e-1;
    end
end

phi(:,:,:,1) = phi_x_interm;
phi(:,:,:,2) = phi_y_interm;
phi(:,:,:,3) = phi_z_interm;

I_m_phi = interpn(X, Y, Z, I_m_orig, phi_x_interm, phi_y_interm, phi_z_interm,'makima'); 

[~,~,~, phi_x_inv,phi_y_inv,phi_z_inv] = PJDC_on_given_mesh3D_LM_short(Mx,My,Mz,phi_x_interm, phi_y_interm, phi_z_interm,X,Y,Z);
toc;

phi_inv(:,:,:,1) = phi_x_inv;
phi_inv(:,:,:,2) = phi_y_inv;
phi_inv(:,:,:,3) = phi_z_inv;

I_f_phi_inv = interpn(X, Y, Z, I_f, phi_x_inv, phi_y_inv, phi_z_inv, 'makima');
% save(['myspace_',name_is,'.mat']);
end
