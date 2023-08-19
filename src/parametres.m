%Autheur : Abel DIDOUH
%Autheur : Gianni PASSANTE

% Fichier d'initialisation de la simulation
% du moteur et de l'asservissement de position par RETOUR D'ETAT 
clear all
close all
Te=0.001;          % Période d'échantillonnage
%********************************************************************************
% Initialisation des parametres du Moteur
%********************************************************************************
R = 5.8;	    % valeurs constructeur
k = 0.0337;
Kv = 0.0497;	% capteur de vitesse V/rd/s
Kp = 10/6.28;	% capteur de position en V/rd
ro=1/20;	    % réducteur mécamique
Jt=5.8743e-04;  % inertie globale

% paramètres du bloc "Frottements secs et visqueux du modele moteur"

Csec =  0.005
Pente = 8.7060e-10

%*************************************************************************************
%	 Gains et Constante de temps
%**************************************************************************************
 [Amod, Bmod, Cmod, Dmod] = linmod('Modele_Lin_BO');
 TM= 3        ;            % constante de temps identifiée  
 GM= 1.5*1/Kv ;            % gain en vitesse identifié        
 Tr= 3*TM/4   ;            % Temps de réponse identifié      
%*************************************************************************************
%	 MATRICES D'ETAT DU MODELE LINEAIRE SANS ACTION INTEGRALE
%**************************************************************************************
Wc = [Bmod, Amod*Bmod];
detWc = det(Wc);

%*************************************************************************************
%	 COMMANDE EN POSITION SANS ACTION INTEGRALE
%**************************************************************************************
% Cahier de charge

 w0= 3/Tr;                 
 zeta= 0.7;
 p1 = -zeta*w0+j*w0*sqrt(1-zeta^2); %Poles complexes
 p2 = -zeta*w0-j*w0*sqrt(1-zeta^2); %Poles complexes
 P = [p1 p2];    

 % Correcteur
 Kc = place(Amod,Bmod,P)
 Ky = Kc*inv(Cmod);


%*************************************************************************************
%	 MATRICES D'ETAT DU MODELE LINEAIRE AVEC ACTION INTEGRALE
%**************************************************************************************

Ai = [Amod, [0;0],
      1     0 0] 

Bi = [Bmod; 0]

Ci = [Cmod, [0;0],
      0 0 ro*Kp]

Di = zeros;


%*************************************************************************************
%	 COMMANDE EN POSITION AVEC ACTION INTEGRALE
%**************************************************************************************

% Cahier de charge
p3 = 5 * (-zeta * w0);
Pi = [p1 p2 p3];       

% Correcteur
Kci = place(Ai,Bi,Pi)
Kyi = Kci*inv(Ci);

%*************************************************************************************
%	 DISCRETISATION DU CORRECTEUR
%**************************************************************************************
% Matrices d'etat du modele discretise
%Te=0.1;
%Te=0.5;
Te = 2;
[Aid, Bid] = c2d(Ai, Bi, Te);
Cid = Ci;
Did = Di;
% Cahier de charge du correcteur discret
% Correcteur discret
Pd = exp(Pi * Te);
Kcid = place(Aid,Bid,Pd);
Kyid = Kcid*inv(Cid);




