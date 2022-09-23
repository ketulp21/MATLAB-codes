E=30e06; % E
nu=0.25; % nu
Ke=zeros(8);
X=[3 5 5 3]; % global x coords
Y=[2 2 4 4]; % global y coords
weight1=[1 1];
weight2=[1 1];
r=[-1/sqrt(3) -1/sqrt(3); 1/sqrt(3) 1/sqrt(3)];
s=[-1/sqrt(3) 1/sqrt(3) ;-1/sqrt(3) 1/sqrt(3)];

for i = 1:2
    for j=1:2
N1 =((r(i,j) + 1)*(s(i,j) + 1))/4;
N2 =-((r(i,j) - 1)*(s(i,j) + 1))/4;
N3 =((r(i,j) - 1)*(s(i,j) - 1))/4;
N4 =-((r(i,j) + 1)*(s(i,j) - 1))/4;

N1r =(s(i,j)/4) + (1/4);
N1s =(r(i,j)/4) + (1/4);
N2r =-(s(i,j)/4) - (1/4);
N2s =(1/4) - (r(i,j)/4);
N3r =s(i,j)/4 - 1/4;
N3s =r(i,j)/4 - 1/4;
N4r =1/4 - s(i,j)/4;
N4s =-(r(i,j)/4) + (-1/4);

J=[N1r*X(1)+N2r*X(2)+N3r*X(3)+N4r*X(4) N1r*Y(1)+N2r*Y(2)+N3r*Y(3)+N4r*Y(4);
  N1s*X(1)+N2s*X(2)+N3s*X(3)+N4s*X(4),N1s*Y(1)+N2s*Y(2)+N3s*Y(3)+N4s*Y(4)];

L=inv(J);

B=[L(1,1)*N1r+L(1,2)*N1s 0 L(1,1)*N2r+L(1,2)*N2s 0 L(1,1)*N3r+L(1,2)*N3s 0 L(1,1)*N4r+L(1,2)*N4s 0;
    0 L(2,1)*N1r+L(2,2)*N1s 0 L(2,1)*N2r+L(2,2)*N2s 0 L(2,1)*N3r+L(2,2)*N3s 0 L(2,1)*N4r+L(2,2)*N4s;
    L(2,1)*N1r+L(2,2)*N1s L(1,1)*N1r+L(1,2)*N1s L(2,1)*N2r+L(2,2)*N2s  L(1,1)*N2r+L(1,2)*N2s L(2,1)*N3r+L(2,2)*N3s L(1,1)*N3r+L(1,2)*N3s L(2,1)*N4r+L(2,2)*N4s L(1,1)*N4r+L(1,2)*N4s ];


D = (E/(1-nu*nu))*[1 nu 0 ; nu 1 0 ;0 0 (1-nu)/2];  %Elasticity Matrix

Ke=B'*D*B*det(J)*weight1(i)*weight2(j)+Ke
    end
end% gauss loop ends