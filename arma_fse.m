function [y,x,w,u] = arma_fse (flag,dim,N,ma,ar,snr,ma_n,ar_n,useed,n1,n2)

% =======================================================================
% function [y,x,w,u] = arma_fse (flag,dim,N,ma,ar,snr,ma_n,ar_n,useed,n1,n2)
% -----------------------------------------------------------------------
% Generates an ARMA sequence (Nx1) corrupted by additive noise
% at SNR=snr. If snr=0, then NO noise is added.
%
% INPUTS 
% dim   = the QAM dimension (4,16,64,...) 
% N     = number of samples desired
% ma    = MA-part of the ARMA channel that generates the process
% ar    = AR-part of the ARMA channel that generates the process
% snr   = Signal-to-noise ratio at the output of the filter 
% ma_n  = MA-part of the filter that filters the noise
% ar_n  = AR-part of the filter for noise generation
% useed = seed value for the random process (for Matlab) (2 X 1)
% n1    = seed value for real part of noise
% n2    = seed value for imaginary part of noise (if it exists)
%
% OUTPUTS
% y     = the output process
% x     = the noise-free process (y - noise)
% w     = the noise process
% u     = the input random process (input to filter)
%
% Author: H. Pozidis,   September 23, 1998
% =======================================================================

if (flag == 'i')
%  disp('Impulse');
  u=[1 zeros(1,N-1)];     % IMPULSE : for testing purposes...
elseif (flag == 'p')
%  disp('PAM-dim');
  u=PAM(dim,N,1,useed(1));
elseif (flag == 'q')
%  disp('QAM-dim');
  u=qam(dim,N,1,useed(1),useed(2));
end

uu=zeros(1,2*N);
uu(1:2:2*N)=u;
u=uu;  clear uu;
x=filter(ma,ar,u);
sx=sum(abs(x).^2)/length(x);
x=x/sqrt(sx);

if (snr > 0)
  if (abs(sum(imag(u)))>0)
    randn('seed',n1);  wr=randn(1,2*N); 
    randn('seed',n2);  wi=randn(1,2*N);  w=wr+i*wi;
    w=filter(ma_n,ar_n,w);  w=w-mean(w);  wr=real(w); wi=imag(w);
    enx=sum(abs(x).^2);  sw=enx/(10^(snr/10));
    enr=sum(wr.^2);  eni=sum(wi.^2);
    wr=wr/(enr^0.5);   wi=wi/(eni^0.5);
    wr=wr*((sw/2)^0.5);  wi=wi*((sw/2)^0.5); w=wr+i*wi;
    y=x+w;
  else
    randn('seed',n1);  w=randn(1,2*N); 
    w=filter(ma_n,ar_n,w);  w=w-mean(w);
    enx=sum(abs(x).^2);  sw=enx/(10^(snr/10));
    enw=sum(w.^2);  
    w=w/(enw^0.5);  
    w=w*(sw^0.5);
    y=x+w;
  end
else
  y=x;  w=0;
end
