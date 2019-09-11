%%%%%%Amandeep Kaur%%%%%
% Assignment-1
% NPM_DICRETE Model 
% Date: 8/2/19

clc; clear variables; close all;
%%%%%%%%%%%%%%%%%%
% Set parameters %
%%%%%%%%%%%%%%%%%%

global alpha beta mu gamma n0 max

alpha = 0.0154;%0.01; % rate of transalation of monomers
beta = 12;%0.9 %rate of conversion of monomers by prion aggregates 
mu = 0.00077; % degregation rate
gamma = 0.00008; % rate of aggregate fragmentation
n0 = 5; %the minimum stable aggregate size

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Initial Conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%time step
t0 = 0; % initial time or left bpound of the interval
tf = 100; % ending time or right bound of the interval
tspan = [t0 tf];%the time interval starting at t0 and ending at tf

%initial Soluble;
x_init = 10;
y_init = 5;

%maximum number of aggregate size
max = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions for the moments and  discrete system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = zeros(1,max)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inital conditions for the Moments %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xMomentInitial = x_init;
% yMomentInitial = 0;
% zMomentInitial = 0;

y0(1)=x_init;

%y0(4) will be the monomer protein
%y0(4)=x_init;

%y(i) is the aggregate concentration of size i for i>=(n0+1)

y0((n0):max)= y_init;

for i=n0:max
    y0(2) = y0(2) + y0(i);
    y0(3) = y0(3) + (i)*y0(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the Discrete System %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tDiscrete,Soln]= ode23(@npm,tspan,y0);%solution of the system

%%%%%%%%%%%%%%%%
% Plot Results %
%%%%%%%%%%%%%%%%

%Compare the moments in the disrete model to the moments from
% the moment equation

calculatedY = zeros(size(tDiscrete));
calculatedZ = zeros(size(tDiscrete));



for i = 1:size(tDiscrete)
    i ;
    for j = n0:max
        j;
        calculatedY(i) = calculatedY(i) + Soln(i,j);
        calculatedZ(i) = calculatedZ(i) + (j)*Soln(i,j);
    end
end



subplot(3,1,1)
%semilogx(tDiscrete,Soln(:,4),'b*-')%DISCRETE SOLUBLE PROTEIN
%hold on
semilogx(tDiscrete,Soln(:,1),'r*-')%MOMENT 
subplot(3,1,2)
semilogx(tDiscrete,calculatedY,'b*-')
hold on
semilogx(tDiscrete,Soln(:,2),'r*-')
subplot(3,1,3)
semilogx(tDiscrete,calculatedZ,'b*-')
hold on
semilogx(tDiscrete,Soln(:,3),'r*-')
