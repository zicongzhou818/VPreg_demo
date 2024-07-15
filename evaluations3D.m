function evaluations3D(I_m, I_f, I_m_phi, I_f_phi_inv, phi, phi_inv,name_is)
    phi_x = phi(:,:,:,1);
    phi_y = phi(:,:,:,2);
    phi_z = phi(:,:,:,3);
    
    phi_x_inv = phi_inv(:,:,:,1);
    phi_y_inv = phi_inv(:,:,:,2);
    phi_z_inv = phi_inv(:,:,:,3);

    mi_fwd = mutualInfo(I_m_phi, I_f);
    mi_bwd = mutualInfo(I_m, I_f_phi_inv);
    mi = mutualInfo(I_m, I_f);
    mi_fwd_r = (mi_fwd-mi)/mi;
    mi_bwd_r = (mi_bwd-mi)/mi;
    thresh = graythresh(I_m_phi);
    BW_I_m_phi = imbinarize(I_m_phi,thresh);

    thresh = graythresh(I_f);
    BW_I_f = imbinarize(I_f,thresh);
    
    j = jaccard(BW_I_m_phi,BW_I_f);
    d = dice(BW_I_m_phi,BW_I_f);
    
    I_m_max=max(I_m,[],"all");
    I_f_max=max(I_f,[],"all");


    ssd_original = sum(((I_m/I_m_max)*I_f_max - I_f).^2,"all");
    r_fwd = sum(((I_m_phi/I_m_max)*I_f_max - I_f).^2,"all")/ssd_original;

    f_m=compute_JD3D(phi_x,phi_y, phi_z);
    jdmax=max(f_m,[],"all");
    jdmin=min(f_m,[],"all");
    
    filename = ['sample2tpl_evaluate_',name_is,'.txt'];
    fileID = fopen(filename, 'w');
    fprintf(fileID,[' Fwd JSC:', num2str(j),' Dice:', num2str(d), ' SSD_r:', num2str(r_fwd), ' MI_incr pct:', num2str(mi_fwd_r*100),'%d\n']);
    fprintf(fileID,[' Phi Jmin:', num2str(jdmin),' Jmax', num2str(jdmax),'%d\n']);

    thresh = graythresh(I_f_phi_inv);
    BW_I_f_psi = imbinarize(I_f_phi_inv,thresh);

    thresh = graythresh(I_m);
    BW_I_m = imbinarize(I_m,thresh);
    
    j2 = jaccard(BW_I_f_psi,BW_I_m);
    d2 = dice(BW_I_f_psi,BW_I_m);
    
    r_bwd = sum(((I_m/I_m_max)*I_f_max - I_f_phi_inv).^2,"all")/ssd_original;

    f_f=compute_JD3D(phi_x_inv, phi_y_inv, phi_z_inv);
    jdmax2=max(f_f,[],"all");
    jdmin2=min(f_f,[],"all");
    

    fprintf(fileID,[' Bwd JSC:', num2str(j2),' Dice:', num2str(d2), ' SSD_r:', num2str(r_bwd), ' MI_incr pct:', num2str(mi_bwd_r*100),'%d\n']);
    fprintf(fileID,[' Phi_inv Jmin:', num2str(jdmin2),' Jmax=', num2str(jdmax2), '%d\n']);
    fclose(fileID);
    disp(['Evauation written to ' filename]);
end