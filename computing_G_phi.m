function [Gx,Gy,Gz]=computing_G_phi(phi1_comp, phi2_comp, phi3_comp, phi_m1, phi_m2, phi_m3, phi_o1,phi_o2,phi_o3, target_x,target_y,target_z)
h=1; tsd=0.001;
[phi_o_jd, ~, ~, ~]=compute_JD_and_Curl3D(phi_o1,phi_o2,phi_o3,h);
%% find P and Q on the current solution phi_comp
Ix = (phi1_comp-target_x);
Iy = (phi2_comp-target_y);
Iz = (phi3_comp-target_z);

[phi_m1_y,phi_m1_x,phi_m1_z]=gradient(phi_m1,h);
[phi_m2_y,phi_m2_x,phi_m2_z]=gradient(phi_m2,h);
[phi_m3_y,phi_m3_x,phi_m3_z]=gradient(phi_m3,h);
[phi_o1_y,phi_o1_x,phi_o1_z]=gradient(phi_o1,h);
[phi_o2_y,phi_o2_x,phi_o2_z]=gradient(phi_o2,h);
[phi_o3_y,phi_o3_x,phi_o3_z]=gradient(phi_o3,h);

[phi_comp_jd, phi_comp_cv1, phi_comp_cv2, phi_comp_cv3]=compute_JD_and_Curl3D(phi1_comp, phi2_comp, phi3_comp,h);

lambda1 = phi_comp_jd;
lambda2x = phi_comp_cv1;
lambda2y = phi_comp_cv2;
lambda2z = phi_comp_cv3;

iix_iiix1 = lambda1.*phi_o_jd.*(phi_m2_y.*phi_m3_z - phi_m2_z.*phi_m3_y) + lambda2y.* phi_o1_z - lambda2z.*phi_o1_y;
iix_iiix2 = lambda1.*phi_o_jd.*(phi_m3_x.*phi_m2_z - phi_m3_z.*phi_m2_x) + lambda2y.* phi_o2_z - lambda2z.*phi_o2_y;
iix_iiix3 = lambda1.*phi_o_jd.*(phi_m2_x.*phi_m3_y - phi_m2_y.*phi_m3_x) + lambda2y.* phi_o3_z - lambda2z.*phi_o3_y;
[~,iix_iiix1_x,~] = gradient(iix_iiix1,h);
[iix_iiix2_y,~,~] = gradient(iix_iiix2,h);
[~,~,iix_iiix3_z] = gradient(iix_iiix3,h);
IIx_IIIx = iix_iiix1_x + iix_iiix2_y + iix_iiix3_z;

iiy_iiiy1 = lambda1.*phi_o_jd.*(phi_m1_z.*phi_m3_y - phi_m3_z.*phi_m2_y) + lambda2z.* phi_o1_x - lambda2x.*phi_o1_z;
iiy_iiiy2 = lambda1.*phi_o_jd.*(phi_m1_x.*phi_m3_z - phi_m1_z.*phi_m2_z) + lambda2z.* phi_o2_x - lambda2x.*phi_o2_z;
iiy_iiiy3 = lambda1.*phi_o_jd.*(phi_m1_y.*phi_m3_x - phi_m1_x.*phi_m3_y) + lambda2z.* phi_o3_x - lambda2x.*phi_o3_z;
[~,iiy_iiiy1_x,~] = gradient(iiy_iiiy1,h);
[iiy_iiiy2_y,~,~] = gradient(iiy_iiiy2,h);
[~,~,iiy_iiiy3_z] = gradient(iiy_iiiy3,h);
IIy_IIIy = iiy_iiiy1_x + iiy_iiiy2_y + iiy_iiiy3_z;

iiz_iiiz1 = lambda1.*phi_o_jd.*(phi_m1_y.*phi_m2_z - phi_m1_z.*phi_m2_y) + lambda2x.* phi_o1_y - lambda2y.*phi_o1_x;
iiz_iiiz2 = lambda1.*phi_o_jd.*(phi_m1_z.*phi_m2_x - phi_m1_x.*phi_m2_z) + lambda2x.* phi_o2_y - lambda2y.*phi_o2_x;
iiz_iiiz3 = lambda1.*phi_o_jd.*(phi_m1_x.*phi_m2_y - phi_m1_y.*phi_m2_x) + lambda2x.* phi_o3_y - lambda2y.*phi_o3_x;
[~,iiz_iiiz1_x,~] = gradient(iiz_iiiz1,h);
[iiz_iiiz2_y,~,~] = gradient(iiz_iiiz2,h);
[~,~,iiz_iiiz3_z] = gradient(iiz_iiiz3,h);
IIz_IIIz = iiz_iiiz1_x + iiz_iiiz2_y + iiz_iiiz3_z;

Gx = Ix + IIx_IIIx*tsd;
Gy = Iy + IIy_IIIy*tsd; 
Gz = Iz + IIz_IIIz*tsd; 
end
