function [Eeff] = Eeff_AAEM(t,t1,q1,q2,q3,q4)
%Eeff_AAEM Age adjusted effective modulus

%The main steps are presented in Bazant and Jirasek (2018): Creep and 
%hygrothermal effects in concrete structures. Equation numbers refer to 
%this book. 

deltat=1.;
deltats = 0.01;
tm = (t+t1)/2.;

%Eq. (4.59)
c1 = 0.08 + 0.0119*log(t1);

%Calculate relaxation function using approximation (4.58)
R = 1/J_B3(t,t1,q1,q2,q3,q4)*(1.+c1*J_B3(t,t1,q1,q2,q3,q4)/...
    (10. * J_B3(t,t-deltat,q1,q2,q3,q4)) * ...
    (J_B3(tm,t1,q1,q2,q3,q4) / J_B3(t,tm,q1,q2,q3,q4)-1))^(-10);

%Calculate t1star in (4.61)
t1star = 0.;
if t1<=t && t<t1+10 * deltats
    t1star = 0.9*t1+0.1*t;
elseif t1 + 10*deltats <t
    t1star = t1+deltats;
end

%Calculate effective modulus Eq. (4.60)
Eeff = (1-R*J_B3(t1star,t1,q1,q2,q3,q4))/(J_B3(t,t1,q1,q2,q3,q4)...
    -J_B3(t1star,t1,q1,q2,q3,q4));

end

