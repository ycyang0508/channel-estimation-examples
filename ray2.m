function h=ray2(rol,delay,mul,support)
% ==========================================================
% function H=RAY2(rol,delay,mul,support);
% ----------------------------------------------------------
%
% Creates a channel as a superposition of 2 scaled and
% delayed raised-cosine pulses. The channel has the form:
%
% h(t) = [a(t,b) + mul*a(t-delay,b)]*w(t);
%
% where a(t,b) is a raised-cosine pulse with roll-off
% factor b (0<b<1), and mul is a multiplicative constant.
% w(t) is a rectangular window of length "support", that
% truncates the length of the channel to "support" samples.
% This is called a "two ray" multipath channel.
%
% Author: H. Pozidis,   September 23, 1998
% ==========================================================

h=rcosine(rol,0)+mul*rcosine(rol,delay);
n=[1:5:201];   %  OVERSAMPLE BY 2
h=h(n);
k=[21-support:21+support-1];
h=h(k);
