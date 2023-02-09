function [J] = J_B3(t,t1,q1,q2,q3,q4)
%Eeff_AAEM Age adjusted effective modulus
%The main equations are in Bazant and Jirasek (2018): Creep and
%Hygrothermal effects in concrete structures.

%Constants page 710
m = 0.5;
n = 0.1;

%Eq. C.3
r = 1.7*t1^0.12 + 9.;

%Eq. C.4
Z = t1^(-m)*log(1+ (t-t1)^n);

%Eq. C.5
Qf = (0.086*t1^(2./9.) + 1.21*t1^(4./9.))^(-1);

%Eq. C2
Q = Qf*(1.+(Qf/Z)^r)^(-1/r);

%Eq. C1
J = q1 + q2*Q + q3*log(1.+(t-t1)^n)+q4*log(t/t1);

end