function [ P_roots ] = MS_root( Y, r, K, N_sectors, NPS, M1, make_plot )
%Mark Wagner, Feb 25, 2020
%Extended root-MUSIC using Manifold separation technique
%Y          = MxL array measurements
%r          = Mx1 sensor positions (line array)
%K          = # DOAs to find
%N_sectors  = # sectors to divide array into
%NPS        = dimension 2 of interpolated manifold about sector
%M1         = dimension 1 of interpolated manifold

[M,~]       = size(Y);                      %dimensions
Ryy         = Y*Y';                         %sample covariance matrix
Ryy(isnan(Ryy) | isinf(Ryy) ) = (1/M^2);
[U,~,~]     = svd(Ryy);                     %svd
En          = U(:,K+1:end);                 %noise subspace
G           = En*En';                       %Grammian'

sectors     = linspace(-1,1,N_sectors);
d_sector    = sectors(2)-sectors(1);
r_ULA       = (0:M1-1)';                    

Proots  = zeros(2*N_sectors*M1,1);
count   = 1;
for i = 2:N_sectors         %Find roots for each sector
    thetas      = linspace(sectors(i)-(.75*d_sector), sectors(i)+(.75*d_sector),NPS);
    A_manifold  = exp(1i*pi.*thetas.*r(:));
    B_manifold  = exp(1i*pi.*thetas.*r_ULA(:));
    G1          = A_manifold*pinv(B_manifold);
    G2          = (G1'*conj(G)*G1);
    poly_coeff  = zeros(2*M1-1,1);
    for j = -M1+1:M1-1
        poly_coeff(j+M1) = sum(diag(G2,j));
    end
    poly_roots  = roots(poly_coeff);   
   
    angle_max   = angle(exp(-1i*pi*(sectors(i)-d_sector/2)));
    angle_min   = angle(exp(-1i*pi*(sectors(i)+d_sector/2)));
    
    if make_plot
        points = 0:.01:2;
        figure(4)
        clf
        zplane(poly_roots)
        hold on
        scatter(real(points.*exp(1i*angle_max)),imag(points.*exp(1i*angle_max)),30,'r','.')
        scatter(real(points.*exp(1i*angle_min)),imag(points.*exp(1i*angle_min)),30,'r','.')
        title('Interpolated array roots')
        axis([-1.5,1.5,-1.5,1.5])
        pause(.1)
    end

    if i < N_sectors    %regular case
        poly_roots(angle(poly_roots)>angle_max | angle(poly_roots)<angle_min | abs(poly_roots) > 1 ) = 0;
        Proots(count:count+length(poly_roots)-1) = poly_roots;
        count = count+length(poly_roots);
    else                %wrap around case
        poly_roots( (angle(poly_roots)>angle_max & angle(poly_roots)<angle_min) | abs(poly_roots) > 1 ) = 0;
        Proots(count:count+length(poly_roots)-1) = poly_roots;
    end
end

P_roots = (sort(Proots(Proots~= 0),'descend'));
P_roots = P_roots(1:K);

if make_plot
    figure(4)
    zplane(P_roots(1:K))
    axis([-1.5,1.5,-1.5,1.5])
    pause(.01)
end

end

