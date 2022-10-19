function [dW1] = Passivity_test_CT(q,dq,B,C,Kd,Kp,e,tau)
%Funzione che permette di calcolare l'energia e la potenza per il calcolo
%della passività del robot
            
        W1 = 0.5*dq'*B*dq;
    
        dW1 = dq'*C*dq - dq'*B*Kd*dq + dq'*B*Kp*e+dq'*tau';
end
