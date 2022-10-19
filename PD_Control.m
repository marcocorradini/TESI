function [PD_Control_results] = PD_Control(q,qA,param,dt,Kp,Kd,maxindex,vertice_inferiore,vertice_superiore,k)
    % Questa funzione funziona come la funzione Computed_torque, con
    % l'aggiunta di una condizione, ovvero che, se l'end effector urta la
    % parete, essa risponde con una forza contraria comportandosi come un
    % corpo elastico (con k elevata)

    %   estraggo i valori dei moventi e i parametri
    m1=param(1,1);
    m2=param(1,2);
    a1=param(1,3);
    a2=param(1,4);
    l1=param(1,5);
    l2=param(1,6);
    I1=param(1,7);
    I2=param(1,8);
    
    % alloco spazio per le variabili usate di seguito 
    E = [];                     % Errore
    dE = [];                    % Derivata errore
    ddQ = zeros(2,1);           % Accelerazione
    dQ = zeros(2,1);            % Velocità
    Q = q;                      % Angoli
    U = zeros(2,1);             % Ingresso controllo 
    TAU = zeros(2,1);           % Coppia applicata dalle forze esterne
    ERR1=zeros(1,1);            % Errore per il calcolo della passività "non compensato" (energia senza 1/2 e' Kp e)
    ERR2=zeros(1,1);            % Errore per il calcolo della passività "compensato" (energia con 1/2 e' Kp e)
    POSITION=zeros(2,1);        % Salva la posizione dell'end effector ad ogni iterata
    
    loop_condition=1;           %condizione necessaria per restare nel loop, se cambia di valore, significa che il robot si ferma ed è arrivato alla posizione finale
    
    e = qA-q;                   % errore iniziale
    de = zeros(2,1);            % derivo l'errore
    dq = zeros(2,1);            % velocia' iniziale
    E = [E, e];                 % vettore dinamico dell'errore
    dE = [dE, de];              % vettore dinamico della derivata dell'errore

    ii=0;
    
    while ii<maxindex &&  loop_condition==1         %condizioni: ii minore del numero massimo di iterazioni consentite e soluzione non ancora trovata
        
        ii=ii+1;
        [xA,yA] = get_EndEffectorPosition(q,param); %estraggo la posizione dell'end effector in base agli angoli theta correnti
        
        if xA>vertice_inferiore(1,1) && yA<vertice_superiore(1,2)        %condizione di interazione con la parete
            dx=[vertice_inferiore(1,1)-xA 0];       %spostamento eleastico della parete
            fext=k*dx;                              %forza di reazione della parete

            tau=fext*[-a2*sin(q(1,1)+q(2,1))-a1*sin(q(1,1)) -a2*sin(q(1,1)+q(2,1));     %calcolo della coppia trasmessa ai giunti
                        a2*cos(q(1,1)+q(2,1))+a1*cos(q(1,1)) a2*cos(q(1,1)+q(2,1))];    %fext*matrice coefficenti, calcolo fatto a mano
        else
            tau=[0 0];                              %se la condizione non è vera, l'end effector non interagiosce con la parete e di conseguenza non si genera alcuna forza
        end

            [B C g]=get_dynamics(q,dq,param);       %estraggo le matrici di inerzia, di Christoffel e di gravità
            u=Kp*e+Kd*de+g;                         %coppia in uscita dal regolatore
               
            % Calcolo l'energia e la potenza in ogni istante per verificare
            % che il controllore sia passivo
            [W1,W2,dW1,dW2] = Passivity_test_PD(q,dq,B,Kd,Kp,e,tau);
    
            err1=dq'*tau'-dW1;                      % Errore per il calcolo della passività "senza compensazione"
            err2=dq'*tau'-dW2;                      % Errore per il calcolo della passività "con compensazione"


            ddq=pinv(B)*(u+tau'-C*dq-g);            %accelerazione calcolata con il modello del robot
            dq=dq+ddq*dt;                           %velocità calcolata integrando l'accelerazione
            q=q+dq*dt;                              %posizione calcolata integrando la velocità

            position=[xA;yA];
            
            e=qA-q;                                 %calcolo l'errore fra la configurazione raggiunta e quella da raggiungere
            de=0-dq;                                %derivo l'errore per il calcolo della coppia del regolatore PD
            
            % Salvo i valori attuali nei vettori dinamici
            E = [E, e];
            dE = [dE, de];
            Q = [Q, q];
            dQ = [dQ, dq];
            ddQ = [ddQ, ddq];
            U = [U, u];
            TAU=[TAU,tau'];
            ERR1=[ERR1,err1];
            ERR2=[ERR2,err2];
            POSITION=[POSITION,position];

        if norm(e)<0.01 && norm(de)<0.001           %errore <j e la velocità finale<m
            loop_condition=0;                       %se le condizioni sono soddisfatte, allora il robot è arrivato ad una configurazione che può andare bene, quindi aggiorno la variabile che fa fermale il loop
        end             
    end
    
    % Salvo i vettori dinamici nella lista corrispondente
    PD_Control_results.E   = E;  
    PD_Control_results.dE  = dE;
    PD_Control_results.ddQ = ddQ;
    PD_Control_results.dQ  = dQ;
    PD_Control_results.Q   = Q;
    PD_Control_results.U   = U;
    PD_Control_results.Tau = TAU;
    PD_Control_results.ERR1= ERR1;
    PD_Control_results.ERR2= ERR2; 
    PD_Control_results.POSITION= POSITION;
end