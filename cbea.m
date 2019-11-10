function [Rec,h]=cbea(snr,LEN,CC,Qq,x1,x2,h1,h2,N1,N2);

% ======================================================================
% function [Rec,h]=cbea(snr,LEN,CC,Qq,x1,x2,h1,h2,N1,N2);
% ----------------------------------------------------------------------
% Implements the CBEA algorithm.
% Reference: H. Pozidis and A.P. Petropulu,
% "Cross-Spectrum Based Blind Channel Identification",
% IEEE Trans. on Signal Processing, vol. 45, no. 12, pp. 2977-2993
%
% INPUTS:
% snr         : signal-to-noise ratio for additive noise processes
% LEN         : number of output symbols to use in the estimation
% CC          : length of the cepstrum window (rectangular here)
% Qq          : length of the autocorrelation window (rectangular here)
% x1,x2       : the observation sequences (single-input, double-output) 
% h1,h2       : the actual sub-channels (for allignement purposes only)
% N1,N2       : (over)-estimates of the lengths of h_min(n), h_max(n)
%
% OUTPUTS:
% Rec         : A matrix containing the estimated h1(n) as its first
%               row, and the estimated h2(n) as its second row.
% h           : the actual channel h(n), from which h1(n) and h2(n) can
%               be obtained by oversampling by a factor of 2. 
%
% Author: H. Pozidis,   September 23, 1998
% ======================================================================

fixlen=128;  L=fixlen;     % the DFT size (preferably a power of 2)
LL=LEN;                   % length of input sequence s(n)
R1=[]; R2=[];    
Q=Qq; B=Qq;   
  
u1=x1(1:LEN);  u2=x2(1:LEN);   % Using only LEN symbols for estimation
Lx1=length(u1); Lx2=length(u2);

T=xcorr(u1,u2,B); [o,p]=max(abs(T));  T=T/T(p);  T=T(:).'; LT=length(T);  
Ry=[T(B+1:LT) zeros(1,fixlen-LT) T(1:B)];
Ry=[Ry(1) fliplr(Ry(2:fixlen))];
Pxx=fft(Ry);               % The cross-spectrum of x1 and x2

R1=xcorr(u1,Q); R1=R1(:).';    
R2=xcorr(u2,Q); R2=R2(:).';
Lr1=length(R1); Lr2=length(R2);
[o,p1]=max(abs(R1)); R1=R1./R1(p1);
[o,p2]=max(abs(R2)); R2=R2./R2(p2);
T1=[R1(p1:Lr1) zeros(1,fixlen-Lr1) R1(1:p1-1)]; 
T1=[T1(1) T1(fixlen:-1:2)];  P1=fft(T1);
T2=[R2(p2:Lr2) zeros(1,fixlen-Lr2) R2(1:p2-1)];
T2=[T2(1) T2(fixlen:-1:2)];  P2=fft(T2);

rcep1=0.5*(ifft(log(abs(P1))));    % REAL CEPSTRUM (Oppenheim-Schafer)
rcep2=0.5*(ifft(log(abs(P2))));
lmin=[1 2*ones(1,fixlen/2-1)];
lmin=[lmin zeros(1,fixlen-length(lmin))];
cep1=rcep1.*lmin; cep2=rcep2.*lmin;
cep1=[cep1(1:CC) zeros(1,fixlen-CC)];
cep2=[cep2(1:CC) zeros(1,fixlen-CC)];
eq1=exp(fft(cep1,fixlen)); eq2=exp(fft(cep2,fixlen));   % The MINIMUM-PHASE
                                % EQUIVALENT SEQUENCES  (Oppenheim-Schafer)

D1=Pxx .* eq1 .* conj(eq2);   
D1=D1(:); pd1=phase(D1);  
psi1=(1/2)*pd1;    psi1=psi1-psi1(1);  % Estimating the phase of hmin(n)

D2=Pxx .* conj(eq1) .* eq2;    
D2=D2(:); pd2=phase(D2);  
psi2=(1/2)*pd2; psi2=psi2-psi2(1);     % Estimating the phase of hmax(n)

%N1=7;                         % An overestimate of the length of hmin(n)
[w1,lser1]=rec_complex(psi1,fixlen,N1);
w1=conj(w1);

%N2=5;                         % An overestimate of the length of hmax(n)
[w2,lser2]=rec_complex(psi2,fixlen,N2);
w2=conj(w2);

rmin = w1;  rmax=w2;
Rmin = roots(rmin); Rmin = Rmin(:).';
Rmax = roots(rmax); Rmax = Rmax(:).';
r1 = poly([Rmin(abs(Rmin)<1) Rmax(abs(Rmax)>1)]);
r2 = poly([1./(Rmin(abs(Rmin)>1)) 1./(Rmax(abs(Rmax)<1))]);
[o,p]=max(abs(r1));  r1=r1/r1(p);      % Normalize w.r.t. max value
[o,p]=max(abs(r2));  r2=r2/r2(p);      % Normalize w.r.t. max value
Rec1=r1;  Rec2=r2;

% ###################################################################
% THE FOLLOWING ARE INCLUDED ONLY FOR ALLIGNEMENT OF THE ESTIMATED
% CHANNELS WITH RESPECT TO THE ACTUAL CHANNEL, FOR PUSPOSES OF
% COMPARISON WHEN PERFORMING MONTE CARLO SIMULATIONS.  THIS IS BECAUSE
% THE CHANNELS ARE RECONSTRUCTED WITHIN AN UNKNOWN CIRCULAR SHIFT.
% ###################################################################

% -------- NORMALIZATION with respect to the MAX value ---------
[o1,p1]=max(abs(h1));  [o2,p2]=max(abs(h2));
Rec1=Rec1*h1(p1);   Rec2=Rec2*h2(p2);
 
%=============================================================%
lr1=length(Rec1); [o,pos1]=max(abs(Rec1));
lh1=length(h1);   [o,ph1]=max(abs(h1));   
zp=10;                                  % Amount of zero-padding

% ----- Set a value for zp, then check if it is big enough ------

MAXNC1=ph1-1+zp;  MAXC1=lh1-ph1+zp;
if ((pos1-1)<MAXNC1) & ((lr1-pos1)<MAXC1)
  Rec1=[zeros(1,MAXNC1-pos1+1) Rec1 zeros(1,MAXC1-lr1+pos1)];
else
  Rec1=zeros(1,2*zp+lh1);  
end
lr2=length(Rec2); [o,pos2]=max(abs(Rec2));
lh2=length(h2);   [o,ph2]=max(abs(h2));
MAXNC2=ph2-1+zp;  MAXC2=lh2-ph2+zp;
if ((pos2-1)<MAXNC2) & ((lr2-pos2)<MAXC2)
  Rec2=[zeros(1,MAXNC2-pos2+1) Rec2 zeros(1,MAXC2-lr2+pos2)];
else
  Rec2=zeros(1,2*zp+lh2);
end
%=============================================================%
Rec=[Rec1;Rec2]; 

