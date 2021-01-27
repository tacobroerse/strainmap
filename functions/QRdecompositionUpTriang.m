function  [U,R,a,b,gamma,theta,Uinv]=QRdecompositionUpTriang(F)
% do upper triangular QR decomposition of deformation gradient into :
% a distortion V and rotation matrix R, via F=R*U
% see Freed  & Srinivasa 2015 (Logarithmic strain and its material
% derivative for a QR decompisition of the deformation gradient)
% notation follows Freed et al. 2020 Laplace Stretch Eulerian and
% Lagrangian formulations

U=zeros(2,2);
R=zeros(2,2);
% right cauchy green deformation tensor
C = F'*F;

%  QT decomposition
% distortion U
U(1,1) = sqrt(C(1,1));
U(1,2) = C(1,2)/U(1,1);
U(2,2) = sqrt(C(2,2)-U(1,2)^2);

% get extensions
% this assumes an order of deformation
% first shear, then extensions
% i.e. F = R * Lambda * gamma
% Lambda = [a 0 ; 0 b] 
% gamma = [1 gamma ; 0 1]
a=U(1,1);
b=U(2,2);
gamma=U(1,2)/U(1,1);

% rotation matrix R = F * inv(U)
Uinv=[1/U(1,1) -U(1,2)/(U(1,1)*U(2,2)) ;...
      0          1/U(2,2) ];
R=F*Uinv;


% rotation angle
theta = atan2(R(2,1),R(1,1));

end

