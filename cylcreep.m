clear;
digits(15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Input parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Units: length in mm, force in N, t in days

%Geometry
ri = 8;
ro = 58;
mr = 100; %steps for r

%Material parameters
E28 = 30000; %target Young's modulus at 28 days for deltat = 0.01days
nu=0.2; %Poisson coefficient
fc28 = 30; %Compressive strength at 28 days
ft28 = 3; %Tensile strength at 28 days

%Parameter used to calculate time dependent tensile strength. 
fib_s = 0.25;

nc = 4.; %Number of cracks
GF = 0.15; %Fracture energy for one crack. Assumed to stay constant
wf = nc*GF/ft28; %%Threshold Exponential stress crack opening curve. For n cracks.

%Corrosion
alpha=2; %corrosion expansion co-ef
icor = 10; %corrosion current density (uA/cm^2) 

t1 = 28; %age of the concrete at the beginning of loading

%Calibration of basic creep parameters according to B3 model
c = 360; %cement content kg/m^3
w = 180; %water content kg/m^3
a = 1860; %aggregate content kg/m^3

%q parameters according to table C.2 in 1/MPa
q1Start = 126.77* fc28^(-0.5);
q2Start = 185.4*c^0.5*fc28^(-0.9);
q3Start = 0.29*(w/c)^4*q2Start;
q4Start = 20.3*(a/c)^(-0.7);

%Multiply values with 1.e-6;
q1Start = q1Start*1.e-6;
q2Start = q2Start*1.e-6;
q3Start = q3Start*1.e-6;
q4Start = q4Start*1.e-6;

%Calculate effective modulus for deltat = 0.01 days at t1 = 28 days
[EeffTest] = Eeff_AAEM( 28.01,28.,q1Start,q2Start,q3Start,q4Start);

%Adjust q parameters
q1 = q1Start*EeffTest/E28;
q2 = q2Start*EeffTest/E28;
q3 = q3Start*EeffTest/E28;
q4 = q4Start*EeffTest/E28;

%Check effective modulus for deltat = 0.01 days at t1 = 28 days
%It should be now equal to E28
[EeffTest] = Eeff_AAEM( 28.01,28.,q1,q2,q3,q4);
if abs(EeffTest-E28)>1.e-6
    error('Something went wrong with calibration of q-values');
end

%Time
mt = 100; % steps for time


%Calculate the end time for prescribed uimax
uimax=5e-2;
delta_t = uimax/(0.0315e-3*(alpha-1)*icor);
t2 = t1+delta_t;

t = linspace(t1+delta_t/mt,t2,mt);


%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creep problem %%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Age-Adjusted Effective Modulus according to AAEM method
Eeff1(mt)= 0.;
for i=1:mt
Eeff1(i) = Eeff_AAEM(t(i),t1, q1,q2,q3,q4); %age-adjusted effective modulus for each t between t(1) and t(mt)
end

ft(mt) = 0.; 
for i=1:mt
ft(i) = fct_MC(t(i),ft28,fc28,fib_s);
end

figure(1);
plot(t,Eeff1/E28,t,ft/ft28);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Determining stress and strain in the cylinder %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initialization of vectors %%%%
r(mr)=0.;
x(mr)=0.;
u_i(mt)=0.;
eps(2,mt,mr)=0.;
sig(2,mt,mr)=0.;
ethcr(mt,mr)=0.;
dethcr(mt,mr)=0.;
eps_tc(2,mr)=0.;
i_tc(mr)=0.;

eps0(2,mt,mr)=0.; %strain without creep
sig0(2,mt,mr)=0.; %stress without creep
ethcr0(mt,mr)=0.;

%Parameters used for the different tests
cond(mt,mr)=0.;
R(10)=0.;
flag(mt,mr)=0.; %flag contains an interger that represents the state of the cylinder :
%1 : uncracked
%2 : cracked

for i=1:mt %time loop
    
	u_i(i) = 0.0315e-3*(alpha-1.)*icor*(t(i)-t1);% Boundary condition imposed at the inner border of the cylinder by the steel expansion
    %steel expansion simplified with q=0 (no transport of the rust into the
    %concrete)

    
    %%%%% Solving the problem with numerical integration of the ODE using bvp4 %%%%%
    
    x = linspace(ri,ro,mr);
    
    solinit = bvpinit(x,@cylinit);
    options = bvpset('RelTol',1.e-10,'AbsTol',1.e-10,'NMax',4*mr); %taken from C.Fahy's code
try    
    sol = bvp4c(@cylode,@cylcreepbc,solinit,options,nu,u_i(i),ri,ro,ft(i),Eeff1(i),wf);
catch ME
    break;
end
    y = deval(sol,x);
%Use solution to compute stresses
    for j=1:mr
        cLength = 2.*pi()*x(j);
        ef = wf/cLength;
        ethcr(i,j) = ...
            crackingstrain(y(1,j),y(2,j),x(j),nu,ft(i),Eeff1(i),ef);

        %Elastic strains
        epsR(i,j) = y(2,j);
        epsT(i,j) = y(1,j)/x(j);
        sigR(i,j) = Eeff1(i)/(1.-nu^2)*(epsR(i,j)+...
            nu*(epsT(i,j)-ethcr(i,j)));
        sigT(i,j) = Eeff1(i)/(1-nu^2)*(epsT(i,j)-ethcr(i,j)+nu*epsR(i,j));
    end

end

%Displaying the internal radial stress as a function of the internal displacement 
sigR_ri = sigR(:,1);
sigR_ri_Norm = sigR_ri*(ri/(ft28*(ro-ri)));

data(:,1) = u_i(:);
data(:,2) = sigR_ri_Norm(:);


save('pressure.dat','data','-ascii')

figure(6);
hold on;
plot(u_i,sigR_ri_Norm,'b');
xlabel('Internal displacement ui');
ylabel('Normalised internal radial stress sig_ri');

