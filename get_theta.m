function [theta1P,theta2P,theta1A,theta2A,O,C1,C2,P1,A1] = get_theta(P,A,param)
    %Restituisce i valori degli angoli iniziali e finali date le coordinate del
    %punto da raggiungere
    
    xP=P(1,1);          %coordinate punto di partenza
    yP=P(2,1);
    xA=A(1,1);          %coordinate punto di arrivo
    yA=A(2,1);
    
    %   estraggo i valori dei moventi e i parametri
    m1=param(1,1);
    m2=param(1,2);
    a1=param(1,3);
    a2=param(1,4);
    l1=param(1,5);
    l2=param(1,6);
    I1=param(1,7);
    I2=param(1,8);
    
    
    theta1P = atan2( (a1^3*yP^2+a1*yP^2*(-a2^2+xP^2+yP^2)-xP*sqrt((-a1^2)*yP^2*(a1^4+(-a2^2+xP^2+yP^2)^2-2*a1^2*(a2^2+xP^2+yP^2))))/(a1^2*yP*(xP^2+yP^2)),(1/(a1^2*(xP^2+yP^2)))*(a1^3*xP+a1*xP*(-a2^2+xP^2+yP^2)+sqrt((-a1^2)*yP^2*(a1^4+(-a2^2+xP^2+yP^2)^2-2*a1^2*(a2^2+xP^2+yP^2)))) );
    
    theta2P = atan2( sqrt((-a1^2)*yP^2*(a1^4+(-a2^2+xP^2+yP^2)^2-2*a1^2*(a2^2+xP^2+yP^2)))/(a1^2*a2*yP),(-a1^2-a2^2+xP^2+yP^2)/(a1*a2) );
                
    theta1A = atan2( (a1^3*yA^2+a1*yA^2*(-a2^2+xA^2+yA^2)-xA*sqrt((-a1^2)*yA^2*(a1^4+(-a2^2+xA^2+yA^2)^2-2*a1^2*(a2^2+xA^2+yA^2))))/(a1^2*yA*(xA^2+yA^2)),(1/(a1^2*(xA^2+yA^2)))*(a1^3*xA+a1*xA*(-a2^2+xA^2+yA^2)+sqrt((-a1^2)*yA^2*(a1^4+(-a2^2+xA^2+yA^2)^2-2*a1^2*(a2^2+xA^2+yA^2)))) );
    
    theta2A = atan2( sqrt((-a1^2)*yA^2*(a1^4+(-a2^2+xA^2+yA^2)^2-2*a1^2*(a2^2+xA^2+yA^2)))/(a1^2*a2*yA),(-a1^2-a2^2+xA^2+yA^2)/(a1*a2) );
    
    
    O=[0;0];
    C1=[a1*cos(theta1P);a1*sin(theta1P)];                                                   % giunto centrale nella configurazione A
    C2=[a1*cos(theta1A);a1*sin(theta1A)];                                                   % giunto centrale nella configurazione B
    P1=[a1*cos(theta1P)+a2*cos(theta1P+theta2P);a1*sin(theta1P)+a2*sin(theta1P+theta2P)];   % punto P per verifica
    A1=[a1*cos(theta1A)+a2*cos(theta1A+theta2A);a1*sin(theta1A)+a2*sin(theta1A+theta2A)];   % punto A per verifica
end