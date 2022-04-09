function [x,y,u,kesi] = dynsim(G,K,deltaT,dx,dy,du,T,x0)
% Simulate a discrete-time plant G with controller K
% deltaT: discrete-time interval
% dx, dy, du: disturbances on x, y, and u, respectively
% x0, system initial state

A = G.A; B = G.B; C = G.C; % D = G.D; assume strictly plant dynamics

nx = size(A,1); % state dimension 
ny = size(C,1); % output dimension
nu = size(B,2); % input dimension
nk = size(K.A,1); % controller state dimension

Tn = floor(T/deltaT);   % number of iterations

x = zeros(nx,Tn+1);
y = zeros(ny,Tn+1);
u = zeros(nu,Tn+1);
kesi = zeros(nk,Tn+1);

x(:,1)    = x0;           % system initial state
kesi(:,1) = zeros(nk,1);  % controller initial state

for i = 1:Tn
    
    % measurement 
    y(:,i)      = C*x(:,i) + dy(:,i);
    
    % controller
    kesi(:,i+1) = K.A*kesi(:,i) + K.B*y(:,i);
    u(:,i)      = K.C*kesi(:,i) + K.D*y(:,i) + du(:,i);
    
    % next state 
    x(:,i+1)    = A*x(:,i) + B*u(:,i) + dx(:,i);
end



end

