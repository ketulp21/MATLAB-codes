function [w,r]=guass_quadrature(n)
if n==2
    w=[1,1;1,1];
    r=[-1/sqrt(3),1/sqrt(3),1/sqrt(3),-1/sqrt(3);
   -1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3)];
elseif n==3
        w=[5/9,8/9,5/9;
        5/9,8/9,5/9];
        r=[-0.77459,0,0.77459,0.77459,0,-0.77459,-0.77459,0,0.77459;
        -0.77459,-0.77459,-0.77459,0,0,0,0.77459,0.77459,0.77459];
end
end