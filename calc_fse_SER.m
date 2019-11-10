function [uest,mmse,symb,err,SER] = calc_fse_SER (fl,dim,real_h,hrec,snr,Lf,ni)

% ==============================================================================
% function [uest,mmse,symb,err,SER] = calc_fse_SER(fl,dim,real_h,hrec,snr,Lf,ni)
%
% dim  =  the QAM dimension (4,16,64,...) (fl = 'q') 
%         or PAM dimension (2,4,8,...)    (fl = 'p')
%
% Author: H. Pozidis,   September 23, 1998
% ==============================================================================

L = 1000;  d = sqrt(dim);
load see; load gaus;
[y,x,w,u]=arma_fse(fl,dim,L,real_h,1,snr,1,1,see(ni,:),gaus(ni,1),gaus(ni,2));

[mmse,mse,f_opt,delta] = MMSE_fse (real_h,hrec,snr,Lf);

uest = conv(y,f_opt);
uest = uest(1:2:2*L);

if (delta <= Lf)
  uest = uest(delta:length(uest));
else
  uest = uest(delta:length(uest));
end
sumsq = sum(abs(uest).^2)/length(uest);
uest = uest/sqrt(sumsq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fl == 'p')
  alphabet=-(dim-1):2:(dim-1);                % values for PAM letters
  cons=alphabet/sqrt(mean(abs(alphabet).^2));
elseif (fl == 'q')
  alphabet = -(d-1):2:(d-1); 	              % values for PAM letters
  [X,Y] = meshgrid(alphabet,alphabet);
  cons = X + Y*j;
  cons = cons(:);
  cons = cons/sqrt(mean(abs(cons).^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ue = uest(:)*ones(1,length(cons));
CONS = ones(length(ue),1)*cons(:).';
eucl = abs(ue-CONS);
[dist,alpha] = min(eucl');

symb = cons(alpha);
u = u(1:2:2*length(uest));
err = u(:) - symb(:);

SER = nnz(err)/length(err);
disp([mmse SER]);

%subplot(121);plot(y,'*');axis('square');
%subplot(122);plot(uest,'*');axis('square');
