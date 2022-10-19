function [dW1] = Passivity_test_CT2(q,dq,B,C,Kd,Kp,e,tau,z,dz,T,alpha)
%Funzione che permette di calcolare l'energia e la potenza per il calcolo
%della passività del robot
            
        W1 = 0.5*dq'*B*dq+T;
    
        dW1 = dq'*C*dq - dq'*B*Kd*dq + alpha*dq'*B*Kp*e + dq'*tau' + z*dz;
end

