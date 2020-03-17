function [ theta_est ] = do_SBLML3( Y, K, r, res )

M           = (size(Y,1));
L           = (size(Y,2));
thetas      = -90:res:90;
r           = r(:);
A           = exp(1i*pi.*sind(thetas).*r);
sigc        = 0.1;           %10.^(-snr/20)*norm(y );
gamma       = ones(length(thetas),1);
gpeak       = 1*ones(K,1);
Ryy         = (Y*Y')./L;
sumRyy      = trace(Ryy);
for j1=1:1000 %:100
    gammaOld        = gamma ;
    GammaMat        = diag(gamma);
    AGA             = A*GammaMat*A';
    SigmaY          = sigc*eye(M)+ AGA;                  %below Wipf (16)
    SigmaX          = GammaMat-GammaMat*A'*inv(SigmaY)*A*GammaMat;
    ApSigmaYinv     = A'/SigmaY;
    mu              = GammaMat*ApSigmaYinv*Y;    %  Wipf {\cal M} (16)
    Sum_mu          = sum(abs(mu.^2),2)/L;
    %  gamma=real(sqrt(Sum_mu./ real(diag(ApSigmaYinv*A))));
 %   gamma=real(sqrt(Sum_mu./ real(diag(A'*inv(SCM)*A))));
 %   gamma=real(Sum_mu./ (1-diag(SigmaX)./diag(GammaMat))); %Wipf (19)
    gamma           = real(Sum_mu +diag(SigmaX)); %Wipf (18)
    %%
    [gpeak,Ilocs]   = findpeaks(gamma,'SORTSTR','descend','NPEAKS',K);
    Am              =A(:,Ilocs);
    Am_inv          =pinv(Am);
    sigc            =((1/L)*sum(sum(abs((Y-A*mu).^2)))+sigc*sum(1-diag(real(SigmaX))./gamma))/M;


    
%     errornorm(j1,isim) = norm(gamma-gammaOld,1)/norm(gamma,1);
  
    if norm(gamma-gammaOld,1)/norm(gamma,1) < 10^(-3)  
        break;
    end
% errorsigma(j1,isim)=sigc;
end
%plot((gamma)), title([num2str(j1) ' ' num2str(sigc) ' ' num2str(errornorm(j1)) ]),

theta_est           =thetas(Ilocs);
% active_set_SBL      =Ilocs;
% Niteration(isim,isnr)=j1;
% sigma(isim,isnr)=sigc;
% GammaEst(isim,isnr,:)=gamma(Ilocs);
% SBLerror(isim,isnr)=errornorm(j1);

end
