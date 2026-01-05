% -----------------------------------------------------------------------
% MODELO DE RAMSEY CON 3 TIPOS DE CHOQUES
% Alfredo Villca
% alfredovillca569@gmail.com
% -----------------------------------------------------------------------

var 
    y       % Producto
    c       % Consumo 
    i       % Inversión 
    k       % Capital físico
    lambda  % Utilidad Marginal del Consumo
    A       % Tecnología
    g       % Gasto público    
    eta     % Preferencias
    ;


varexo 
    e       % Choque de tecnología
    eg      % Choque de gasto público
    eeta    % Choque de preferencias
    ;

parameters 
    lambda_ss 
    c_ss 
    k_ss 
    y_ss 
    i_ss 
    alpha 
    beta 
    delta 
    sigma 
    rho 
    rhog 
    rhoeta 
    sigmae 
    sigmaeg 
    sigmaeeta
    ;

alpha     = 0.36;
beta      = 0.98;
delta     = 0.025;

sigma      = 1.5;
rho        = 0.95;
rhog       = 0.9;
rhoeta     = 0.9;
sigmae     = 0.01;
sigmaeg    = 0.01;
sigmaeeta  = 0.01;


A_ss      = 1;
g_ss      = 0;
eta_ss    = 1;
k_ss      = ((alpha*beta)/(1 - beta*(1-delta)) )^(1/(1-alpha)); 
y_ss      = k_ss^alpha;
i_ss      = delta*k_ss;
c_ss      = y_ss - (i_ss + g_ss);
lambda_ss = c_ss^(-sigma) ;




model;
[name = 'Utilidad marginal del consumo']
eta*c^(-sigma) = lambda;

[name = 'Ecuacion del Euler']
lambda = alpha*beta*lambda(+1)*A(+1)*k^(alpha-1) + beta*(1-delta)*lambda(+1); 

[name = 'Demanda agregada']
y = c + i + g;

[name = 'Funcion de produccion']
y = A*k(-1)^alpha;

[name = 'Ley de movimiento del capital']
k = i + (1-delta)*k(-1);

[name = 'Tecnologia']
log(A) = rho*log(A(-1)) + e;

[name = 'Gasto publico']
g = rhog*g(-1)+eg;

[name = 'Preferencias']
log(eta) = rhoeta*log(eta(-1)) + eeta;

end;

initval;
lambda	=   lambda_ss;
c	      =	  c_ss;
k	      =	  k_ss;
y       =   y_ss;
i       =   i_ss;
A       =   A_ss;
g       =   g_ss;
eta     =   eta_ss;
end;

steady;
resid;

shocks;
var e;    stderr (sigmae);
var eg;   stderr (sigmaeg);
var eeta; stderr (sigmaeeta);
end;

stoch_simul(order=1, irf=100, periods=500); 

