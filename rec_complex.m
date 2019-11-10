function [rec,lse] = rec_complex (theta,NF,N)

% ===================================================================
% function [rec,lse] = rec_complex (theta,NF,N)
% -------------------------------------------------------------------
% Reconstructs a sequence based on knowledge of its phase within a
% linear phase component. The length of the sequence is unknown.
%
% theta   :  the phase sequence
% NF      :  the FFT size
% N       :  an estimate of the length of the reconstructed sequence
%
% rec     :  the reconstructed sequence
% lse     :  least-squares-error
%
% Author: H. Pozidis,   September 23, 1998
% ===================================================================

lse=[];  NI=N;  W=-1;   % Range within which to search from the length N.
vv=[NI-W:NI+W]; vv=4*vv-3;  vv=sum(vv);
REC=zeros(vv,NI+abs(W));
nloop=1;

for N=NI-W:NI+W
  for shift=0:0.5:0.5
    psi1=theta+(2*pi/NF)*[0:NF-1]*shift;
    psi1=psi1(:);

    for a=0:N-1
      vec1=[-(N-1-a):a];   
      pies1=(2*pi/NF)*[0:NF-1]'*vec1;
      for j=1:length(vec1)
        ang1(1:NF,j)=psi1; 
      end
      Phi1=sin(pies1+ang1);
      Phi2=-cos(pies1+ang1);  
      PHI=[Phi1 Phi2];

      [U,S,V]=svd(PHI); 
      [si,sj]=sort(diag(S));
      UU=V(:,sj);
      u=UU(:,1);

      rec=u(1:N)+sqrt(-1)*u(N+1:2*N);
      rec=rec(:).';
      err=(norm(PHI*u)^2)/NF;
      lse=[lse err];

      REC(nloop,:)=[rec zeros(1,NI+W-length(rec))];
      nloop=nloop+1;
    end
  end
end

[lsm,lsi]=min(lse); rec=REC(lsi,:);
f=find(rec ~= 0); rec=rec(f); %disp(size(rec));
