function [ Tu,count ] = AP_Gridless( Y,r,K,max_iter,tol, make_plot )
%Mark Wagner, Dec 31 2019
%extended Alternating Projection based gridless beamforming (e-APG) for
%arbitrary array geometry. Given measurements Y, and number of sources K,
%estimate DOAs. max_iter is maximum number of iterations, tolerance is a
%break condition calling the function to end angle between subspaces of Tu
%and the projection of its extended Vandermonde matrix is low enough.

% M = number of sensors
% L = number of snapshots
% K = number of sources
% Y = MxL measurement matrix
% r = vector of sensor positions
% Tu= MxM hermitian toeplitz matrix of rank K
% Z = LxL Hermitian symmetric matrix
% max_iter = maximum number of iterations
% tol = break condition

%% end initialization
Y           = (Y./norm(Y,'fro'));
r           = r(:).';
[M,L]       = size(Y);                  %number of sensors, number of snapshots
ULA         = 0;
if diag(squareform(pdist(r')),1) == ones(M-1,1)
    ULA     = 1;
end
theta_ADN   = find_CBF_peaks(Y, r, 1, 1/M, 1e-10);  %find maxima of CBF, (also known as atomic dual norm)
V           = exp(1i*2*pi*theta_ADN(:)'.*r(:));     
Tu          = (1/M)*(V*V');
Z           = zeros(L); 

count           = 0;
Bs              = zeros(max_iter,1);        %allocate space to record break conditions
S_prev          = zeros(M+L);
for i = 1:max_iter
    S                   = [Tu,Y;Y',Z];                      %construct S
    S                   = PSD(S);                           %PSD projection
    Tu                  = S(1:M,1:M);                       %new Tu
    if ULA
        Tu              = Tproj(Tu);                        %Toeplitz projection
        Bs(i)           = norm(S-S_prev,'fro');
        S_prev          = S;
    else
        [Tu,~,~]        = ITP(K,r,Tu);                      %Irregular Toeplitz projection
        Z               = S(M+1:end,M+1:end);
        Bs(i)           = norm(S-S_prev,'fro');
        S_prev          = S;
        if Bs(i) < tol
            count =         i;
            break
        end
    end
    
end

if count == 0
    count = max_iter;
end

if make_plot
    figure(3)
    semilogy(abs((Bs)))
    grid on
    title('$\| S^{(k)}- S^{(k+1)} \|_F$ vs. Iteration', 'Interpreter','latex')
    xlabel('Iteration')
    ylabel('$\| S^{(k)}- S^{(k+1)} \|_F$', 'Interpreter','latex')
    pause(.01)
end
end

%% Revmoved Nesterov acceleration
%     Tu_prev             = Tu;
%
%     if accel
%         k                   = (sqrt(1+4*k_prev^2) + 1 )/2;
%         Tu                  = Tu + ((k_prev-1)/k)*(Tu-Tu_prev);
%         k_prev              = k;
%     end
