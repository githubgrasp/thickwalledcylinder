function [ft] = fct_MC(t, ft28, fc28, fib_s)
%fct_Mc Calculate time dependence of tensile strength based on Model Code
%2010

fcm = exp(fib_s * ( 1. - sqrt(28./t) ) ) * fc28;

%Calculate the aged tensile strength
if fcm >= 58.
    fib_ft = 2.12 * log(1. + 0.1 * fcm);
elseif fcm <= 20.
    fib_ft = 0.07862 * fcm; % 12^(2/3) * 0.3 / 20 = 0.07862
else
   fib_ft = 0.3 * (fcm - 8.)^(2./3.); %5.1-3a
end

%Calculate the 28 day tensile strength
if fc28 >= 58.
    fib_ft28 = 2.12 * log(1. + 0.1 * fc28);
elseif fc28 <= 20.
    fib_ft28 = 0.07862 * fc28; % 12^(2/3) * 0.3 / 20 = 0.07862
else
    fib_ft28 = 0.3 * (fc28 - 8.)^(2./3.); %5.1-3a
end

ft = fib_ft/fib_ft28 * ft28;

end

