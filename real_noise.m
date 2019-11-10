function [x1,x2,n1,n2]=real_noise(x1,x2,snr,s1,s2)

% ===================================================
% function [x1,x2,n1,n2]=real_noise(x1,x2,snr,s1,s2)
%
% Adds noise to processes x1(n) and x2(n), at 
% SNR=snr, using seeds s1 and s2.
%
% Author: H. Pozidis,   September 23, 1998
% ===================================================

if (nargin == 5)
 if (snr < 0)
  n1=0;  n2=0;
 else
  randn('seed',s1);  n1=randn(size(x1));  
  randn('seed',s2);  n2=randn(size(x2));
  n1=n1-mean(n1); n2=n2-mean(n2);

  ss1=sum(abs(x1).^2);  ss2=sum(abs(x2).^2);
  sw1=ss1/(10^(snr/10));sw2=ss2/(10^(snr/10));

  n1=n1/(norm(n1));    n2=n2/(norm(n2));
  n1=n1*(sw1^0.5); n2=n2*(sw2^0.5);
  x1=x1+n1;    x2=x2+n2;
 end
elseif (nargin == 3)
  t1=x2;
  t2=snr;
  snr=t1;
  s1=t2;
  if (snr < 0)
    x2=0;
  else
    randn('seed',s1);  n1=randn(size(x1));
    n1=n1-mean(n1);

    ss1=sum(abs(x1).^2); 
    sw1=ss1/(10^(snr/10));

    n1=n1/(norm(n1));    
    n1=n1*(sw1^0.5); 
    x1=x1+n1;
    x2=n1;
  end
end

clear ss1 ss2 sw1 sw2 s1 s2;
