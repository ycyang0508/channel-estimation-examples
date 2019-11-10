function s = qam (M,N_s,sigma,sd1,sd2)
 
% =====================================================================
% function s = qam (M,N_s,sigma,sd1,sd2)
% ---------------------------------------------------------------------
%
% creates a sequence of N_s samples of M-QAM with variance sigma^2,
% where M is of the form 4*N^2 (for integer N).
%
% ----- Downloaded from the SPIB (Rice University) -----
% =====================================================================

 clear j
 % a convenient quantity
 N = sqrt(M/4);					% M = 4N^2

 % check if valid M (via checking if N is an integer)
 if N ~= floor(N),
   error('M must be of the form 4*N^2 for integer N'); 
 end

 % generate M-QAM with constellation points at +/-{1,3,5,7,etc.}
 rand('seed',sd1);
 sr = 2*(floor(2*N*rand(1,N_s))-N+0.5);
 rand('seed',sd2);
 si = 2*(floor(2*N*rand(1,N_s))-N+0.5);
 s = sr + j*si; 

 % adjust variance 
 var_cur = 2/N*sum((1+2*[0:N-1]).^2);		% current variance
 s = s*(sigma/sqrt(var_cur));
