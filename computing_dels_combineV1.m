function [del_F,del_f, del_g]=computing_dels_combineV1(X_0, I_m_local,I_f,phi_X,N,f)
I_m_temp = interpn(X_0(:,:,:,1),X_0(:,:,:,2),X_0(:,:,:,3),I_m_local,phi_X(:,:,:,1),phi_X(:,:,:,2),phi_X(:,:,:,3),'makima'); 
Img_diff = I_m_temp-I_f; 
[I_m_y, I_m_x, I_m_z] = gradient(I_m_temp); 
lap_B(:,:,:,1) = Img_diff(2:N-1,2:N-1,2:N-1).*I_m_x(2:N-1,2:N-1,2:N-1);
lap_B(:,:,:,2) = Img_diff(2:N-1,2:N-1,2:N-1).*I_m_y(2:N-1,2:N-1,2:N-1);
lap_B(:,:,:,3) = Img_diff(2:N-1,2:N-1,2:N-1).*I_m_z(2:N-1,2:N-1,2:N-1);
parfor ith=1:3
    del_F(:,:,:,ith) = pois3fft(lap_B(:,:,:,ith));
end

[B1_x2,  B1_x1,  B1_x3] = gradient(del_F(:,:,:,1)); 
[B2_x2,  B2_x1,  B2_x3] = gradient(del_F(:,:,:,2)); 
[B3_x2,  B3_x1,  B3_x3] = gradient(del_F(:,:,:,3)); 
del_f = - f.*(B1_x1 + B2_x2 + B3_x3);
del_g(:,:,:,1) = - B3_x2 + B2_x3;
del_g(:,:,:,2) =   B3_x1 - B1_x3;
del_g(:,:,:,3) = - B2_x1 + B1_x2;
end