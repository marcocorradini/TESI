function [B,C,g] = get_dynamics(q,dq,param)

    % Raggruppo le funzioni  get_B,  get_C e  get_g per avere solo un'unica funzione nel main
    
    B = get_B(q,param);
    C = get_C(q,dq,param);
    g = get_g(q,param);
end

