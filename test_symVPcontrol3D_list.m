clc
clear
close all

folderpath = '/home/joe/Documents/SymVPcontrol3D_demo';
ls_all = dir(folderpath);
keys = '_01_MR_T1w_MPR_corrected_bet_lin_on_mni152_Warped.nii.gz';
ls_img = ls_all(contains({ls_all.name},keys));
ls_img_name = {ls_img.name};
tpl_is = strcat(folderpath,'/mni_icbm152_t1_tal_nlin_sym_09c_masked.nii.gz');
I_f = niftiread(tpl_is);
I_f_info = niftiinfo(tpl_is);


J_thr = 0.2;
regrid_bd = 4;
aa=length(ls_img_name);

parfor ii=1:aa
    I_m_is = strcat(folderpath, '/', ls_img_name{ii});
    [~, name_is, ~] = fileparts(I_m_is);
    [~, name_is, ~] = fileparts(name_is);
    I_m = single(niftiread(I_m_is));
    I_m_info = niftiinfo(I_m_is);
    I_m_info.Datatype='single';
    
    [phi,I_m_phi,phi_inv,I_f_phi_inv]=Img_Reg3D_InvConsbyVP_lap_combineV1_parallel(I_m,I_f,J_thr,regrid_bd,name_is);
    evaluations3D(I_m, I_f, I_m_phi, I_f_phi_inv, phi, phi_inv);

    niftiwrite(I_m_phi,[name_is,'_SymVPcont_Mapped'],I_f_info,"Compressed",true);
    niftiwrite(phi,[name_is,'_SymVPcont_Map'],"Compressed",true);
    niftiwrite(I_f_phi_inv,[name_is,'_SymVPcont_InverseMapped'],I_m_info,"Compressed",true);
    niftiwrite(phi_inv,[name_is,'_SymVPcont_InverseMap'],"Compressed",true);
    
end



seg01_wm = single(niftiread('mni_icbm152_wm_tal_nlin_sym_09c.nii.gz'));
seg02_gm = single(niftiread('mni_icbm152_gm_tal_nlin_sym_09c.nii.gz'));
seg03_csf = single(niftiread('mni_icbm152_csf_tal_nlin_sym_09c.nii.gz'));

[Mx,My,Mz] = size(I_f);
[X, Y, Z] = ndgrid(1:Mx,1:My,1:Mz);
phi_inv=niftiread('Sample_01_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMap.nii.gz');
phi_x_inv=phi_inv(:,:,:,1);
phi_y_inv=phi_inv(:,:,:,2);
phi_z_inv=phi_inv(:,:,:,3);

sample01_info = niftiinfo('Sample_01_T1w_corrected_bet_lin_on_mni152_Warped.nii.gz');
sample01_info.Datatype='single';

sample01_seg01_wm_phi_inv = interpn(X, Y, Z, seg01_wm, phi_x_inv, phi_y_inv, phi_z_inv, 'makima');
sample01_seg01_wm_phi_inv=(sample01_seg01_wm_phi_inv-min(sample01_seg01_wm_phi_inv,[],"all"))/max(sample01_seg01_wm_phi_inv,[],"all");
niftiwrite(sample01_seg01_wm_phi_inv,'Sample_01_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMapped_wm',sample01_info,"Compressed",true);

sample01_seg02_gm_phi_inv = interpn(X, Y, Z, seg02_gm, phi_x_inv, phi_y_inv, phi_z_inv, 'makima');
sample01_seg02_gm_phi_inv=(sample01_seg02_gm_phi_inv-min(sample01_seg02_gm_phi_inv,[],"all"))/max(sample01_seg02_gm_phi_inv,[],"all");
niftiwrite(sample01_seg02_gm_phi_inv,'Sample_01_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMapped_gm',sample01_info,"Compressed",true);

sample01_seg03_csf_phi_inv = interpn(X, Y, Z, seg03_csf, phi_x_inv, phi_y_inv, phi_z_inv, 'makima');
sample01_seg03_csf_phi_inv=(sample01_seg03_csf_phi_inv-min(sample01_seg03_csf_phi_inv,[],"all"))/max(sample01_seg03_csf_phi_inv,[],"all");
niftiwrite(sample01_seg03_csf_phi_inv,'Sample_01_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMapped_csf',sample01_info,"Compressed",true);

psi_inv=niftiread('Sample_02_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMap.nii.gz');
psi_x_inv=psi_inv(:,:,:,1);
psi_y_inv=psi_inv(:,:,:,2);
psi_z_inv=psi_inv(:,:,:,3);

sample02_info = niftiinfo('Sample_02_T1w_corrected_bet_lin_on_mni152_Warped.nii.gz');
sample02_info.Datatype='single';

sample02_seg01_wm_psi_inv = interpn(X, Y, Z, seg01_wm, psi_x_inv, psi_y_inv, psi_z_inv, 'makima');
sample02_seg01_wm_psi_inv=(sample02_seg01_wm_psi_inv-min(sample02_seg01_wm_psi_inv,[],"all"))/max(sample02_seg01_wm_psi_inv,[],"all");
niftiwrite(sample02_seg01_wm_psi_inv,'Sample_02_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMapped_wm',sample02_info,"Compressed",true);

sample02_seg02_gm_psi_inv = interpn(X, Y, Z, seg02_gm, psi_x_inv, psi_y_inv, psi_z_inv, 'makima');
sample02_seg02_gm_psi_inv=(sample02_seg02_gm_psi_inv-min(sample02_seg02_gm_psi_inv,[],"all"))/max(sample02_seg02_gm_psi_inv,[],"all");
niftiwrite(sample02_seg02_gm_psi_inv,'Sample_02_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMapped_gm',sample02_info,"Compressed",true);

sample02_seg03_csf_psi_inv = interpn(X, Y, Z, seg03_csf, psi_x_inv, psi_y_inv, psi_z_inv, 'makima');
sample02_seg03_csf_psi_inv=(sample02_seg03_csf_psi_inv-min(sample02_seg03_csf_psi_inv,[],"all"))/max(sample02_seg03_csf_psi_inv,[],"all");
niftiwrite(sample02_seg03_csf_psi_inv,'Sample_02_T1w_corrected_bet_lin_on_mni152_Warped_SymVPcont_InverseMapped_csf',sample02_info,"Compressed",true);





% F-mod comp#:0 ei/ti:100/188 ts:1.296e-06 r/rr:0.72421/0.72421 J_min:0.63145 rMI:%3.4481
% F-mod comp#:0 ei/ti:100/187 ts:1.1782e-06 r/rr:0.74259/0.74259 J_min:0.54655 rMI:%4.4889
%  r_fwd:0.74193 J_min_comp:0.54762
% F-mod comp#:1 ei/ti:25/43 ts:1.3306e-06 r/rr:0.72307/0.97458 J_min:0.94979 rMI:%4.7696
%  r_fwd:0.72421 J_min_comp:0.63145
% F-mod comp#:1 ei/ti:27/47 ts:1.3041e-06 r/rr:0.69715/0.96263 J_min:0.94379 rMI:%3.4946
%  r_fwd:0.7003 J_min_comp:0.62197
% f&g-mod comp#:2 ei/ti:12/12 ts:1.7959e-05 r/rr:0.5466/0.78053 J_min:0.19106 rMI:%4.0523
%  r_fwd:0.72316 J_min_comp:0.51838
% f&g-mod comp#:2 ei/ti:16/16 ts:2.1829e-05 r/rr:0.56269/0.7781 J_min:0.17508 rMI:%5.1913
%  r_fwd:0.56746 J_min_comp:0.14404
% f&g-mod comp#:3 ei/ti:22/30 ts:1.7559e-05 r/rr:0.49303/0.86883 J_min:0.1656 rMI:%4.1033
%  r_fwd:0.55286 J_min_comp:0.19833
% f&g-mod comp#:3 ei/ti:32/57 ts:1.1959e-05 r/rr:0.42902/0.77599 J_min:0.18413 rMI:%4.1191
%  r_fwd:0.47976 J_min_comp:0.032705
% f&g-mod comp#:4 ei/ti:43/79 ts:1.1634e-05 r/rr:0.39858/0.83078 J_min:0.1962 rMI:%6.3025
%  r_fwd:0.43256 J_min_comp:0.099568
% f&g-mod comp#:4 ei/ti:43/79 ts:1.1634e-05 r/rr:0.36411/0.84177 J_min:0.19326 rMI:%4.8524
%  r_fwd:0.40386 J_min_comp:0.0083458
%  r_fwd:0.368 J_min_comp:0.064758
% Phi_inv found ts:0.7778 r:8.9096e-09 ei/ti:82/93
% Elapsed time is 10353.367302 seconds.
%  r_fwd:533.1387 r_bwd:453.9688 regrid:4 elasped:10447.6777
%  Fwd JSC:0.94241 Dice:0.97035 SSD_r:0.40386 MI_incr: %7.1412
%  Phi Jmin:0.0083458 Jmax3.1656
%  Bwd JSC:0.96915 Dice:0.98433 SSD_r:0.38956 MI_incr: %13.7505
%  Phi_inv Jmin:0.20007 Jmax=54.1029
% Phi_inv found ts:1.3991 r:6.0261e-09 ei/ti:76/76
% Elapsed time is 10402.328038 seconds.
%  r_fwd:770.4171 r_bwd:721.3793 regrid:4 elasped:10497.2969
%  Fwd JSC:0.94371 Dice:0.97104 SSD_r:0.368 MI_incr: %6.0455
%  Phi Jmin:0.064758 Jmax2.6069
%  Bwd JSC:0.97223 Dice:0.98592 SSD_r:0.37275 MI_incr: %14.3166
%  Phi_inv Jmin:0.33908 Jmax=12.6338




