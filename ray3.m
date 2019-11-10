function h=ray3(rol,delay,mul,support)
% ===================================================================
% function h=ray3(rol,delay,mul,support)
% -------------------------------------------------------------------
%
% Creates a channel as a superposition of 3 scaled and
% delayed raised-cosine pulses. The channel has the form:
%
% h(t)=[a(t,b)+mul(1)*a(t-delay(1),b)+mul(2)*a(t-delay(2),b)]*w(t);
%
% where a(t,b) is a raised-cosine pulse with roll-off
% factor b (0<b<1), and mul is a multiplicative constant.
% w(t) is a rectangular window of length "support", that
% truncates the length of the channel to "support" samples.
% This is called a "three ray" multipath channel.
%
% Author: H. Pozidis,   September 23, 1998
% ===================================================================

d1=delay(1);  d2=delay(2);
m1=mul(1);    m2=mul(2);

h=rcosine(rol,0)+m1*rcosine(rol,d1)+m2*rcosine(rol,d2);
n=[1:5:201];
k=[21-support:21+support-1];
h=h(n);
h=h(k);
