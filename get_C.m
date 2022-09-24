function [C] = get_C(q,dq,param)
    
    %   estraggo i valori dei moventi e i parametri
    theta1=q(1,1);
    theta2=q(2,1);
    
    dtheta1=dq(1,1);
    dtheta2=dq(2,1);

    m1=param(1,1);
    m2=param(1,2);
    a1=param(1,3);
    a2=param(1,4);
    l1=param(1,5);
    l2=param(1,6);
    I1=param(1,7);
    I2=param(1,8);

    %   date le posizioni dei moventi e i parametri del sistema, trova la matrice C
    C=zeros(2,2);
    C=[-m2*a1*l2*sin(theta2)*dtheta2 -m2*a1*l2*sin(theta2)*(dtheta1+dtheta2);
        m2*a1*l2*sin(theta2)*dtheta1 0];
end