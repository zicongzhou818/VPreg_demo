function [JD_phi, CV_phi1, CV_phi2, CV_phi3]=compute_JD_and_Curl3D(Phi1,Phi2,Phi3,h)
[Phi1y,Phi1x,Phi1z]=gradient(Phi1,h);
[Phi2y,Phi2x,Phi2z]=gradient(Phi2,h);
[Phi3y,Phi3x,Phi3z]=gradient(Phi3,h);


JD_phi=Phi1x.*(Phi2y.*Phi3z-Phi2z.*Phi3y)...
    -Phi2x.*(Phi1y.*Phi3z-Phi1y.*Phi3z)...
    +Phi3x.*(Phi1y.*Phi2z-Phi1z.*Phi2y);

CV_phi1=Phi3y-Phi2z;
CV_phi2=Phi1z-Phi3x;
CV_phi3=Phi2x-Phi1y;
end