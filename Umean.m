% This function calculates the mean velocity and the friction factor
% This function receives:
% deq,S,g,hc,CD,a,phi,Cmean,lambda,alpha,beta,C,kappa,Pi as defined in Run.m
% This function returns:
% Depth-aeraged mean velocity (Umean); friction factor (f)

function [Umean,f] = Umean(deq,S,g,hc,CD,a,phi,Cmean,lambda,alpha,beta,C,kappa,Pi,UUD,y0)
%% Substitute the solved deq to its relevant parameters
Y1 = deq/lambda; %elevation of the location where c=0.01 (m)
Y90 = deq/(1-Cmean); %elevation of the location where c=0.90 (m)
zeta=1.6*deq/hc; %dimensionless inflection point velocity (zeta = ui/uUD)(-)
um=sqrt(g*S*(deq-hc)); %shear velocity at the canopy top (m/s)
utot=sqrt(g*S*(deq+hc*(phi-1))); %shear velocity at the channel bottom (m/s)
UFSY1=sqrt(2*g*S/(CD*a))+(zeta*UUD-UUD)*(1+tanh((Y1-beta*hc)/(alpha*hc)))+um/kappa*log((Y1-beta*hc)/(y0))+um*2*Pi/kappa*sin(pi*Y1/(2*deq/lambda))^2; %constant velocity in the free-surface layer
%% Calculation of the mean velocity of each layer
UUDmean = (deq/(lambda*Y90))*sqrt(2*g*S/(CD*a)); %uniform distribution layer
UMLmean = (deq/(lambda*Y90))*sqrt(2*g*S/(CD*a))*(zeta-1)*(1+(lambda*alpha*hc/deq)*log(cosh(((beta/alpha)-deq/(lambda*alpha*hc)))/cosh(beta/alpha))); %mixing layer
ULLmean = (1/Y90)*(um/kappa)*((deq/lambda)-(beta*hc))*(log((deq/lambda - beta*hc)/hc)+kappa*C-1); %log layer
UWFmean = um*deq/(lambda*Y90)*Pi/kappa; %wake function layer
UFSmean= UFSY1*deq/Y90*(1-1/lambda); %free-surface layer

%% Calculation of the mean velocity
Umean = UUDmean+UMLmean+ULLmean+UWFmean+UFSmean; %summation of all the mean velocity components

%% Calculation of the friction factor
f=8*(utot/Umean)^2; 


end



    