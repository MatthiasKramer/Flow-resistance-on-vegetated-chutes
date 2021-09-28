%% Solver for flow resistance on grassed channels 
%Version 1.0 (September 2021)
%Developed by: Matthias Kramer and Hanwen Cui
%Works with Matlab 2021a
%Contact:
%m.kramer@adfa.edu.au

%When using this code, please cite/refer to the following reference:
%--------------------------------------------------------------------------%
%H. Cui, S. Felder, and M. Kramer (2019)
%Five-layer velocity model predicts flow resistance of aerated flows down grass-lined spillways 
%Journal of Hydraulic Engineering
%--------------------------------------------------------------------------%

clear all
close all

%% Input parameters
q=0.3; %define the flow rate that needs to be calculated
X = ['q (input) = ', num2str(q),' (m^2/s)']; disp (X) 
syms deq %equivalent clear water depth (m)

% Channel characteristics
S=sind(10.8); %chute slope (-)
g=9.81; %(m/s^2)

% Grass properties
hc=0.02; %deflected canopy height (m)
CD=1; %drag coefficient (-)
a=281; %frontal area per volume (1/m)
phi=1; %porosity (-)

% Aeration properties 
Cmean = 0.1; %depth-averaged air concentration (-) 0 for single-phase flows
lambda=1.2; %dimensionless elevation of the location where c=0.01 (lambda=deq/Y1) (-) 1 for single-phase flows
Y1=deq/lambda; %elevation of the location where c=0.01 (m)
Y90=deq/(1-Cmean); %elevation of the location where c=0.90 (m)

% Velocity profile fitting parameters
alpha=0.09; %dimensionless mixing layer length scale (alpha=Le/hc) (-)
beta=0.85; %dimensionless elevation of the inflection point (beta=yi/hc) (-)
C=3.4; %integration constant log-law, 5.2 for subcritical flows (-) 
kappa=0.41; %von Karman constant (-)
Pi=0.3; %wake function term (-)
UUD=sqrt(2*g*S/(CD*a)); %in-canopy uniform velocity;
y0=hc*exp(-kappa*C); %hydraulic roughness(m)
zeta=1.6*deq/hc; %dimensionless inflection point velocity (zeta = ui/uUD)(-)

%% Numerical Solution for deq
um=sqrt(g*S*(deq-hc)); %shear velocity at the canopy top(m/s)
eqnLeft=q/deq; %continuity equation
eqnRight=(deq/(lambda*Y90))*sqrt(2*g*S/(CD*a))+(deq/(lambda*Y90))*sqrt(2*g*S/(CD*a))*(zeta-1)*(1+(lambda*alpha*hc/deq)*log(cosh(((beta/alpha)-deq/(lambda*alpha*hc)))/cosh(beta/alpha)))+(1/Y90)*(um/kappa)*((deq/lambda)-(beta*hc))*(log((deq/lambda - beta*hc)/hc)+kappa*C-1)+um*deq/(lambda*Y90)*Pi/kappa+(sqrt(2*g*S/(CD*a))+(zeta*UUD-UUD)*(1+tanh((Y1-beta*hc)/(alpha*hc)))+um/kappa*log((Y1-beta*hc)/(y0))+um*2*Pi/kappa*sin(pi*Y1/(2*deq/lambda))^2)*deq/Y90*(1-1/lambda); %Eq. (11) in ......

% Numerical solver for deq
Sol = vpasolve(eqnLeft == eqnRight, deq);
deqSol= double(Sol);

% Calculation of the mean velocity and friction factor
[UmSol,fSol]=Umean(deqSol,S,g,hc,CD,a,phi,Cmean,lambda,alpha,beta,C,kappa,Pi,UUD,y0);
deqSolhc=deqSol/hc; %submergence of the interested grass-lined flow

X = ['d (solution) = ', num2str(deqSol),' (m)']; disp (X) 
X = ['u (solution) = ', num2str(UmSol),' (m/s)']; disp (X) 
X = ['f (solution) = ', num2str(fSol), ' (-)']; disp (X) 
