% Load mat variables
data = load("modelProfilesForQgLinStab");
depth= data.depth;
rho  = data.rho;
uVel = data.vVel;
vVel = data.uVel;
bigF = data.bigF;
beta = data.beta;
betaT= [0,0];

% initialize choice of wavenumbers
myVecK = logspace(-3.5,-1.5,1000);
myVecL = logspace(-4.5,-2.5,1000);

%% 

% Call qggrz function
[wiMax, wrMax, psiVec] = qggrz(depth,rho,uVel,vVel,bigF,beta,betaT,myVecK,myVecL,1);

%% 

% test
max(wrMax,[],"all")
max(wiMax,[],"all")
contourf(myVecL,myVecK,wiMax)
set(gca, "XScale", "log")
set(gca, "YScale", "log")