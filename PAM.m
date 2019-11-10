function s = PAM (M,N_s,sigma,sd)

% ==================================================================
% function s = PAM (M,N_s,sigma,sd)
%
% creates a sequence of N_s samples of M-PAM with variance sigma^2,
% where M is an even integer.
%
% ---- Downloaded from the SPIB (Rice University) ----
% ==================================================================


 % check if valid M 
 if M/2 ~= floor(M/2),
   error('M must be an even integer'); 
 end

 % generate M-PAM with constellation points at +/-{1,3,5,7,etc.}
 rand('seed',sd); 
 s = 2*(floor(M*rand(1,N_s))-M/2+0.5); 

 % adjust variance 
 var_cur = 2/M*sum((1+2*[0:M/2-1]).^2);		% current variance
 s = s*(sigma/sqrt(var_cur));
