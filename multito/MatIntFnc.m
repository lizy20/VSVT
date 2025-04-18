function [w,dwdy,V,dVdy] = MatIntFnc(y,xi,beta,type,penal)
switch(type)
  case('SIMP') 
    h=y;
    w = y.^penal;
    V = y;
    dwdy = penal*y.^(penal-1);
    dVdy = ones(size(y));

  case('SIMP-H')
    h =  ( tanh(xi*beta) + tanh(beta*(y-xi))) /(tanh(xi*beta)+ tanh(beta*(1-xi)));     %1-exp(-beta*y)+y*exp(-beta) ;
    w = h.^penal;
    V = h; % h;
    dhdy = beta*(1-(tanh(beta*(y-xi))).^2 ) /(tanh(xi*beta)+tanh(beta*(1-xi))) ;           %beta*exp(-beta*y)+exp(-beta); %
    dwdy = penal*h.^(penal-1).*dhdy;
    dVdy = dhdy;
end
end