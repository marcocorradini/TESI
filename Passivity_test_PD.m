function [W1,W2,dW1,dW2] = Passivity_test_PD(q,dq,B,Kd,Kp,e,tau)
%Funzione che permette di calcolare l'energia e la potenza per il calcolo
%della passività del robot
    
        W1 = 0.5*dq'*B*dq;
        W2= W1 + 0.5*e'*Kp*e;
    
        dW1 = -dq'*Kd*dq - dq'*Kp*e+dq'*tau';
        dW2 = -dq'*Kd*dq+dq'*tau';
end