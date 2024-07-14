function U=pois3fft(F)
%
% Solve 3-d Poisson equation
%
%      -div grad u = f on [0,1]x[0,1]x[0,1]
%
% with 0 boundary condition, i.e., u({0,1},y,z)=u(x,{0,1},z)=u(x,y,{0,1})=0.
%
% Method:
%
%   1) Set up the grid by partitioning along x-, y-,and z-directions
%      equidistant:
%
%         h=1/(n+1); xi=i*h; yj=j*h, zk=k*h
%
%   2) Let U(n-by-n-by-n) and F(n-by-n-by-n) with U(i,j,k)=u(xi,yj,zk) to be determined
%      and F(i,j,k)=f(xi,yj,zk). Then U satisfies
%
%           U *_1 Tn + U *_2 Tn + U *_3 Tn = h^2*F,
%
%      where Tn(n-by-n)=tridiag(-1,2,-1).
% 
% 
%
%   3) Solve the linear system by FFT (direct sine transformation)
%
% Ref: http://www.ec-securehost.com/SIAM/ot56.html, section 6.7
%
%    RCL 11/15/2009
%
% Input
%
%     F   matrix (n-by-n-by-n)
%         F(i,j,k)=f(xi,yj,zk)
%
% Output
%
%     U   matrix (n-by-n-by-n)
%         U(i,j,k) approximates u(xi,yj,zk)

[n,n1,n2]=size(F);
if n1~=n || n2~=n
   disp('pois3ftt: F must be cubic, i.e., n-by-n-by-n');
   return
end
%h=1/(n+1); 
h=1;
G=h^2*F; 
U=zeros(n,n,n);

for k=1:n    
    Uwk=dst(G(:,:,k).',n); Uwk=Uwk.'; U(:,:,k)=dst(Uwk,n);
end
for i=1:n
    % Uwk=U(i,:,:)'
    Uwk=zeros(n,n);
    for j=1:n
        Uwk(:,j)=U(i,j,:);
    end
    Uwk=dst(Uwk,n);
    U(i,:,:)=Uwk.';
end

theta=(pi/(2*(n+1)))*(1:n); tmp=4*(sin(theta')).^2;
% G nolonger needed; used as work array;
Uwk=tmp*ones(1,n)+ones(n,1)*tmp';
for k=1:n
    G(:,:,k)=Uwk+tmp(k);
end
U=U./G;

for k=1:n    
    Uwk=idst(U(:,:,k).',n); Uwk=Uwk.'; U(:,:,k)=idst(Uwk,n);
end
for i=1:n
    % Uwk=U(i,:,:)'
    Uwk=zeros(n,n);
    for j=1:n
        Uwk(:,j)=U(i,j,:);
    end
    Uwk=idst(Uwk,n);
    U(i,:,:)=Uwk.';
end
