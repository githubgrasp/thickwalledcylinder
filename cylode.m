function dydr = cylode(x,y,nu,~,~,~,ft,E,wf)
     cLength = 2.*pi()*x;
     ef = wf/cLength;
     [ethcr] = crackingstrain(y(1),y(2),x,nu,ft,E,ef);
    
     %Use ethcr to compute A 
     paramA = 1.-ft/E*(1-nu^2)/ef*exp(-ethcr/ef);
 
     dydr = [y(2)
           -y(2)/x*(paramA-nu)/(paramA-nu^2) + ...
           y(1)/x^2*(paramA-nu)/(paramA-nu^2) - ...
           ethcr/x*paramA*(1-nu)/(paramA-nu^2)];
end