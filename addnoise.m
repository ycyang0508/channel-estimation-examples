function [rx1,ix1,rx2,ix2]=addnoise(x1,x2,snr,s1,s2,s3,s4)

% =================================================================
% function [rx1,ix1,rx2,ix2]=addnoise(x1,x2,snr,s1,s2,s3,s4)
%
% Generates noise sequences and adds them to the output sequences
%
% x1  :  1st output sequence
% x2  :  2nd output sequence
% snr :  the signal-to-noise-ratio
% s1  :  seed for randn  (also: s2,s3,s4)
%
% rx1 :  real part of noisy output sequence x1
% rx2 :  real part of noisy output sequence x2 
% ix1 :  imaginary part of noisy output sequence x1
% ix1 :  imaginary part of noisy output sequence x2
%
% REMARKS:
% If snr<0, no noise is added --- NOISE-FREE case
% If only 4 input arguments are supplied, then the output sequence
% is assumed real-valued, and real-valued noise is added.
%
% Author: H. Pozidis,   September 23, 1998
% =================================================================

if (nargin == 7)
 if (snr >= 0)
  rx1=real(x1); ix1=imag(x1); rx2=real(x2); ix2=imag(x2);
  randn('seed',s1);  nr1=randn(size(rx1));  
  randn('seed',s2);  ni1=randn(size(ix1));
  randn('seed',s3);  nr2=randn(size(rx2));  
  randn('seed',s4);  ni2=randn(size(ix2));
  nr1=nr1-mean(nr1); ni1=ni1-mean(ni1);
  nr2=nr2-mean(nr2); ni2=ni2-mean(ni2);

  ss1=sum(abs(x1).^2);  ss2=sum(abs(x2).^2);
  sw1=ss1/(10^(snr/10));sw2=ss2/(10^(snr/10));

  sr1=sum(nr1.^2);      si1=sum(ni1.^2);
  sr2=sum(nr2.^2);      si2=sum(ni2.^2);
  nr1=nr1./(sr1^0.5);   ni1=ni1./(si1^0.5);
  nr2=nr2./(sr2^0.5);   ni2=ni2./(si2^0.5);

  nr1=nr1*((sw1/2)^0.5); ni1=ni1*((sw1/2)^0.5);
  nr2=nr2*((sw2/2)^0.5); ni2=ni2*((sw2/2)^0.5);

  rx1=rx1+nr1;    ix1=ix1+ni1;
  rx2=rx2+nr2;    ix2=ix2+ni2;
  clear nr1 ni1 nr2 ni2 ss1 ss2 sw1 sw2 sr1 sr2 si1 si2;
 else
  rx1=real(x1);  ix1=imag(x1);  rx2=real(x2);  ix2=imag(x2);
 end
elseif (nargin == 4)
  t1=x2; t2=snr; t3=s1;
  snr=t1; s1=t2; s2=t3;
  if (snr >= 0)
    rx1=real(x1); ix1=imag(x1); 
    randn('seed',s1);  nr1=randn(size(rx1));
    randn('seed',s2);  ni1=randn(size(ix1));
    nr1=nr1-mean(nr1); ni1=ni1-mean(ni1);
   
    ss1=sum(abs(x1).^2);  sw1=ss1/(10^(snr/10));
	
    sr1=sum(nr1.^2);      si1=sum(ni1.^2);
    nr1=nr1./(sr1^0.5);   ni1=ni1./(si1^0.5);
	 
    nr1=nr1*((sw1/2)^0.5); ni1=ni1*((sw1/2)^0.5);
	  
    rx1=rx1+nr1;    ix1=ix1+ni1;
    clear nr1 ni1 ss1 sw1 sr1 si1;
  else
	rx1=real(x1);  ix1=imag(x1);
  end
else
  error('The routine ''addnoise.m'' takes 4 or 7 arguments');
end
