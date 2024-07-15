clc
clear
close all

folderpath = '/home/joe/Documents/SymVPcontrol3D_demo';
ls_all = dir(folderpath);
keys = '_T1w_corrected_bet_lin_on_mni152_Warped.nii.gz';
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
    I_m_info.Datatype = 'single';
    
    [phi, I_m_phi, phi_inv, I_f_phi_inv] = ImgReg_SymVPcontrol(I_m, I_f, J_thr, regrid_bd, name_is);
    evaluations3D(I_m, I_f, I_m_phi, I_f_phi_inv, phi, phi_inv, name_is);

    niftiwrite(I_m_phi, [name_is, '_SymVPcont_Mapped'], I_f_info, "Compressed", true);
    niftiwrite(phi, [name_is, '_SymVPcont_Map'], "Compressed", true);
    niftiwrite(I_f_phi_inv, [name_is, '_SymVPcont_InverseMapped'], I_m_info, "Compressed", true);
    niftiwrite(phi_inv, [name_is, '_SymVPcont_InverseMap'], "Compressed", true);
end




