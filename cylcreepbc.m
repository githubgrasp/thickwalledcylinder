function res = cylcreepbc(ya,yb,nu,ui,ri,ro,ft,E,wf)

     cLength = 2.*pi()*ro;
     ef = wf/cLength;
     ethcr = crackingstrain(yb(1),yb(2),ro,nu,ft,E,ef);

res = [ ya(1) - ui
        (yb(2) + nu * (yb(1)/ro - ethcr))/(1.-nu^2)];