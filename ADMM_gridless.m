function [Tu] = ADMM_gridless(Y_temp,tau,iter,nesterov)
%NOTE: nesterov PARAMETER TAKES BINARY VALUE 1 OR 0 TO TURN ON OR OFF
%NESTEROV MOMENTUM

%size of mesurements
Y_size = size(Y_temp);
n = Y_size(1,1);
L = Y_size(1,2);

%parameters
rho = 1.0; % penalty parameter in augmented Lagrangian

nesterov_momentum = nesterov;  % set to 1 to use Nesterov schedule for overrelaxation.
momentum = 1.0;         % (constant momentum.  Only used if Nesterov is off.)

%initialization
ZOld = zeros(n+L);                      %actually the S matrix
Lambda = zeros(n+L); 
normalizer = 1./[n; 2*((n-1):-1:1)'];
e1 = zeros(n,1); e1(1)=1;
theta = 1;

for count = 1:iter
    
    %update the variables W,X, and u (X = measurements, u = T(u), W = Z)
    X = ((rho)/(1+2*rho))*(Y_temp/rho...
        + ZOld(1:n,n+1:n+L)/2 + ZOld(n+1:n+L,1:n)'/2 ...
        - Lambda(1:n,n+1:n+L)/2/rho-Lambda(n+1:n+L,1:n)'/2/rho);
    u = normalizer.*(toeplitz_adjoint(ZOld(1:n,1:n) - Lambda(1:n,1:n)/rho)  - (tau/2)*n/rho*e1);
    W = ZOld(n+1:n+L,n+1:n+L)-Lambda(n+1:n+L,n+1:n+L)/rho - (tau/2)/rho*eye(L);
    
    %temp is the matrix that should be psd.
    temp = [toeplitz(u),X; X', W];
    
    % update momentum using nesterov rule
    if nesterov_momentum
        momentum = 1+theta*(1/theta -1);
        theta = (sqrt(theta^4+4*theta^2)-theta^2)/2;
    end
    
    %projection of Q onto the semidefinite cone
    Q = (momentum*temp+Lambda+(1-momentum)*ZOld);
    [V,E] = eig((Q+Q')/2);
    e = diag(E);
    idx = (e>0);
    Z = V(:,idx)*diag(e(idx))*V(:,idx)';
    Z = (Z+Z')/2;
    
    %stop criteria.
%     pri_res = temp - Z;
%     
%     dual_res = rho*[sum(Z(1:n,n+1:n+L)-ZOld(1:n,n+1:n+L),2);
%         toeplitz_adjoint(Z(1:n,1:n)-ZOld(1:n,1:n));
%         sum(Z(n+1:n+L,n+1:n+L)-ZOld(n+1:n+L,n+1:n+L),2)];
%     dual_var_adj = rho*[sum((Lambda(1:n,n+1:n+L)+Lambda(n+1:n+L,1:n)')/2,2);
%         (toeplitz_adjoint(Lambda(1:n,1:n)));
%         sum(Lambda(n+1:n+L,n+1:n+L),2)];
%     
%     pri_tol = sqrt(size(Z,1)*size(Z,2))*tol_abs + tol_rel*max(norm(temp,'fro'),norm(Z,'fro'));
%     dual_tol = sqrt(length(dual_var_adj))*tol_abs + tol_rel*norm(dual_var_adj);
%     
%     err_rec_primal = norm(pri_res,'fro');
%     err_rec_dual = norm(dual_res);
%     
%     converged = and(err_rec_primal<pri_tol,err_rec_dual<dual_tol);
%     if converged, break; end
    
    Lambda = Lambda + momentum*temp+(1-momentum)*ZOld-Z;
    ZOld = Z;
end


Tu = toeplitz(u);

end


function T = toeplitz_adjoint(A)
N = size(A,1);
T = zeros(N,1);
T(1) = sum(diag(A));
for n = 1:(N-1)
    T(n+1) = 2*sum(diag(A,n));
end
end
