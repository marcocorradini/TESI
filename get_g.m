function [g] = get_g(q,param)
    
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
    G=param(1,9);

    %   date le posizioni dei moventi e i parametri del sistema, trova la matrice g
    g=zeros(2,1);
    g=[(m1*l1+m2*a1)*G*cos(theta1)+m2*G*l2*cos(theta1+theta2);
        m2*G*l2*cos(theta1+theta2)];
end