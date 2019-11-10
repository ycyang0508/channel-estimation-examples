function y = rcosine (roloff,tau)

% ======================================================
% function y=rcosine(roloff,tau)
%
% roloff  : the roloff factor for the raised-cosine
%           pulse (a number between 0 and 1)
% tau     : the multipath delay (avoid integer numbers
%           since they create identifiability problems 
%
% Author: H. Pozidis,   September 23, 1998
% ======================================================

T=1;
t=[0.1:0.1:10*T];   Lt=length(t);
tau=tau*T;   delay=ceil(tau/0.1);

x=sin(pi*t/T)./(pi*t/T);
x=x.*(cos(roloff*pi*t/T)./(1-(4*(roloff^2)*(t.^2)/(T^2))));

y=[1 x];
y=[fliplr(x) y];

y=[zeros(1,delay) y(1:2*Lt+1-delay)];
