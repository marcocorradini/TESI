function [B] = get_B(q,param)
    
    %   estraggo i valori dei moventi e i parametri
    theta1=q(1,1);
    theta2=q(2,1);

    m1=param(1,1);
    m2=param(1,2);
    a1=param(1,3);
    a2=param(1,4);
    l1=param(1,5);
    l2=param(1,6);
    I1=param(1,7);
    I2=param(1,8);

    %   date le posizioni dei moventi e i parametri del sistema, trova la matrice B 
    B=zeros(2,2);
    B=[m1*l1^2+I1+m2*(a1^2+l2^2+2*a1*l2*cos(theta2))+I2 m2*(l2^2+a1*l2*cos(theta2))+I2;
        m2*(l2^2+a1*l2*cos(theta2))+I2 m2*l2^2+I2];
end