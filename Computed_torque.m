function [q,y] = Computed_torque(q,qA,param,dt,Kp,Kd,tau,maxindex)
    % Questa funzione integra al suo interno un regolatore PD che permette di generare una
    % coppia U (in base all'errore tra la configurazione attuale e quella a 
    % cui si vuole arrivare) che viene poi utilizzata per calcolare, mediante 
    % l'equazione del modello del robot, l'accelerazione. Integrando l'accelerazione
    % si trovano le velocità e la posizione step per step e inserendo tutto
    % in un loop si arriva alla configurazione finale mediante il metodo di
    % Eulero
    
    % alloco spazio per le variabili usate di seguito 
    e=zeros(2,1);                               %errore
    ddq=zeros(2,1);                             %accelerazione
    dq=zeros(2,1);                              %velocità
    
    loop_condition=1;                           %condizione necessaria per restare nel loop, se cambia di valore, significa che il robot si ferma ed è arrivato alla posizione finale
    
    e=qA-q;                                     %errore iniziale
    de=0-dq;                                    %derivo l'errore
    
    ii=0;
    
    while ii<maxindex &&  loop_condition==1     %condizioni: ii minore del numero massimo di iterazioni consentite e soluzione non ancora trovata
        
        ii=ii+1;
        
        [B C g]=get_dynamics(q,dq,param);       %estraggo le matrici di inerzia, di Christoffel e di gravità
        U=C*dq+g+B*(Kp*e+Kd*de);                %coppia in uscita dal regolatore
               
        ddq=pinv(B)*(U+tau-C*dq-g);             %accelerazione calcolata con il modello del robot
        dq=dq+ddq*dt;                           %velocità calcolata integrando l'accelerazione
        q=q+dq*dt;                              %posizione calcolata integrando la velocità
        
        y(ii,:)=q';                             %salvo in un vettore y le posizioni parziali raggiunte
        
        e=qA-q;                                 %calcolo l'errore fra la configurazione raggiunta e quella da raggiungere
        de=0-dq;                                %derivo l'errore per il calcolo della coppia del regolatore PD
        
        
        if norm(e)<0.01 && norm(de)<0.001       %errore <j e la velocità finale<m
            loop_condition=0;                   %se le condizioni sono soddisfatte, allora il robot è arrivato ad una configurazione che può andare bene, quindi aggiorno la variabile che fa fermale il loop
        end
    end
end

