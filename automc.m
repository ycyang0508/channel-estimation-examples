function [recm,h,recm1,recm2,mmse,ser]=automc(snr,LEN,CC,Qq)

% --------------------------------------------------------------
% function [recm,h,recm1,recm2,mmse,ser]=automc(snr,LEN,CC,Qq)
% --------------------------------------------------------------
%
% Monte Carlo simulations of the CBEA algorithm.
% 
% snr  :  the SNR of the additive noise
% LEN  :  the number of symbols of the output signals
% CC   :  the length of the cepstrum window
% Qq   :  the length of the autocorrelation window
%
% recm :  the vector of reconstructed channels h(n)
% h    :  the actual channel
% recm1:  the vector of reconstructed channels h1(n)
% recm2:  the vector of reconstructed channels h2(n)
%
% Author: H. Pozidis,   September 23, 1998
% --------------------------------------------------------------

NI=50;  dim=2;
Lmin=5;  Lmax=7;   % Example lengths for the intermediate sequences
slen=LEN;  no_seg=LEN/slen;
recm=[]; recm1=[]; recm2=[];

h=ray2(0.11,.7,.8,6);  Lh=length(h); h=reshape(h,2,Lh/2);
%h=ray3(0.12,[0.4 1.2],[0.6 0.8],6);  Lh=length(h); h=reshape(h,2,Lh/2);
h1=h(1:2:Lh);  h2=h(2:2:Lh);

for k=1:NI
  randn; seeds(k,1)=randn('seed'); randn; seeds(k,2)=randn('seed');
  randn; gaus(k,1)=randn('seed');  randn; gaus(k,2)=randn('seed');
  randn; gaus(k,3)=randn('seed');  randn; gaus(k,4)=randn('seed');
end

for jk=1:NI
  s1=seeds(jk,1); s2=seeds(jk,2);
  s3=gaus(jk,1); s4=gaus(jk,2); s5=gaus(jk,3); s6=gaus(jk,4);
  s=PAM(dim,LEN,1,s1);  s=s(:);
%  s=qam(dim,LEN,1,s1,s2);  s=s(:);
  x1=conv(s,h1);  x2=conv(s,h2);
  if (abs(sum(imag(s))) > 0.001)
    [rx1,ix1,rx2,ix2] = addnoise(x1,x2,snr,s3,s4,s5,s6);
    x1=rx1+sqrt(-1)*ix1;  x2=rx2+sqrt(-1)*ix2;
  else
    [x1,x2]=real_noise(x1,x2,snr,s3,s4);
  end

 for kk=1:no_seg
   u1=x1(slen*(kk-1)+1:slen*kk);  u2=x2(slen*(kk-1)+1:slen*kk);
   r=cbea(snr,slen,CC,Qq,u1,u2,h1,h2,Lmin,Lmax);   % for QAM inputs
   HR(kk,:)=r(:).';
 end
 if (no_seg == 1)
   recm(jk,:)=HR;
 else
   recm(jk,:)=mean(HR);
 end 

  [ues,mmse(jk),sym,er,ser(jk)]=calc_fse_SER('p',dim,h(:),recm(jk,21:32),snr,32,jk); 
  recm1(jk,:)=recm(jk,1:2:length(recm(jk,:)));
  recm2(jk,:)=recm(jk,2:2:length(recm(jk,:)));
  disp(jk);
end
