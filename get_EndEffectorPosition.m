function [x1,y1] = get_EndEffectorPosition(q,param)
    % Presi in ingresso i valori dei moventi, calcola le coordinate dell'end effector
    
    %   estraggo i valori dei moventi e i parametri
    m1=param(1,1);
    m2=param(1,2);
    a1=param(1,3);
    a2=param(1,4);
    l1=param(1,5);
    l2=param(1,6);
    I1=param(1,7);
    I2=param(1,8);

    theta1=q(1,1);
    theta2=q(2,1);

    x1=a1*cos(theta1)+a2*cos(theta1+theta2);
    y1=a1*sin(theta1)+a2*sin(theta1+theta2);
end