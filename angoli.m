function [SOLA] = angoli(A,B,param)
    
    xA=A(1,1);
    yA=A(1,2);
    xB=B(1,1);
    yB=B(1,2);
    
    m1=param(1,1);
    m2=param(1,2);
    a1=param(1,3);
    a2=param(1,4);
    l1=param(1,5);
    l2=param(1,6);
    I1=param(1,7);
    I2=param(1,8);

SOLA=solve([a1*cos(ang1A)+a2*cos(ang1A+ang2A)==xA,a1*sin(ang1A)+a2*sin(theta1A+ang2A)==yA],[ang1A,ang2A]);
end

