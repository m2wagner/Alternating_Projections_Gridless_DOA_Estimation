function [Tu_next] = ADMM_NUA(Y,r,K,max_iter, tau)
%NOTE: 
r       = r(:);
[M,L]   = size(Y);
%initialization
theta_ADN   = find_CBF_peaks(Y, r, 1, 1/M, 1e-10);  %find maxima of CBF, (also known as atomic dual norm)
V           = exp(1i*2*pi*theta_ADN(:)'.*r(:));     
Tu          = (1/M^2)*(V*V');
S           = [Tu,Y;Y',0];
Lambda      = zeros(M+L); 

for count = 1:max_iter
    %update Y, Tu, and Z
    Y_next  = (1/3)*(Y +  (S(1:M,M+1:M+L)/2 + S(M+1:M+L,1:M)'/2) - (Lambda(1:M,M+1:M+L)/2+Lambda(M+1:M+L,1:M)'/2) );
    if diag(squareform(pdist(r)),1) == ones(M-1,1)
        Tu_next         = Tproj( S(1:M,1:M) - Lambda(1:M,1:M) );
    else
        [Tu_next,~,~]   = ITP(K,r,S(1:M,1:M)-tau*Lambda(1:M,1:M));
    end
    Z_next  = S(M+1:M+L,M+1:M+L) - tau*(Lambda(M+1:M+L,M+1:M+L));
    
    temp    = [Tu_next,Y_next; Y_next', Z_next];    
    Q       = temp+Lambda;
    Z       = PSD(Q);
    Lambda  = Lambda + temp -Z;
    S       = Z;
end

end



