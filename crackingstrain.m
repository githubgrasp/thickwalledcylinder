function [ethcr] = crackingstrain(y1,y2,x,nu,ft,E,ef)
%Function to determine the cracking strain iteratively from a exponential
%stress-crack opening curve, using a Newton method

if E*(y1/x + nu*y2)/(1.-nu^2.) <= ft
    ethcr = 0; 
else
    ethcr = 0;
    R = 1;
    iter = 0;
    while(abs(R)>1.e-6)
        if(iter >100)
            print('newton method did not converge');
        end

        R = ft*exp(-ethcr/ef) - E*(y1/x - ethcr + nu*y2 )/(1.-nu^2); %Residual which has to be zero        
        dRdethcr = -ft*exp(-ethcr/ef)*1/ef + E/(1-nu^2); %derivative of residual
        ethcr = ethcr - R/dRdethcr; %General newton equation        
        iter=iter+1;
    end
end

end