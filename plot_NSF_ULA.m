%% Run Alternating Projections based Gridless beamforming for ARBITRARY array geometry
clear; clc; 
rng_seed = 0;
rng(rng_seed)

L                       = 10;                       %number of samples
M                       = 21;                       %sensors
K                       = 3;                        %number of signals
SNR                     = 100;                       %dB
SRC_RANGE               = 1;
RANDOM_SENSORS          = 0;

%% Generate necessary arrays
sin_thetas              = sort(mygen_DOAs(K,1/M)*2-1);  %sin(K random DoAs)
sin_thetas              = [-.126,.275,.67];
DoA_rad                 = asin(sin_thetas);           %angles, radians
DoA_deg                 = DoA_rad/(pi)*180;               %angles, degrees
DoA_phasors             = exp(1i*pi*sin_thetas);        %complex angles
if RANDOM_SENSORS
    r                       = [1:M-1]+1*rand(1,M-1);            %generate random sensor spacing
    r                       = [0,r]';                           %first sensor position is zero
else
    r = [0:M-1]';
end

%% Generate Signals
[ Y, Y2, noise, A, A2,X]= gen_signals_SNR( L,M,K,SNR,DoA_phasors,r,SRC_RANGE );
Y2                      = Y2+noise;
Ryy                     = Y2*Y2';
[U,~,~]                 = svd(Ryy);
En                      = U(:,K+1:end);
G                       = En*En';
% fun                     = @(z) (null_spec_polynomial( G, r, exp(1i*pi*sin(angle(z))) )); 
fun                     = @(z) (null_spec_polynomial( G, r, z )); 

%% Evaluate over z plane
xy_min                  = -1.5;
xy_max                  = abs(xy_min);
samples                 = 1000;                              %number of samples
f                       = linspace(-.5,.5,samples);         %coarse grid of frequencies between +/-.5
x                       = xy_min:.01:xy_max;
y                       = xy_min:.01:xy_max;
[X,Y]                   = meshgrid(x,y);
Z                       = X+(1i*Y);
Z_circ                  = exp(1i*2*pi.*f);
S_circ                  = zeros(1,length(Z_circ));
for i = 1:length(Z_circ)
    S_circ(i)                  = abs(fun(Z_circ(i)));
end
S                       = zeros(size(Z)); 
for i = 1:length(x)
    for j = 1:length(y)
        S(i,j)          = 10*log10(abs(fun(Z(i,j))));
    end
end
lsamp                   = 110;
levels                  = linspace(-100,10*log10(4*M),lsamp);
levels                  = 10.^(levels./10);




figure(1)
clf
subplot(1,2,1)
[C,h] =contour(X,Y,S,levels,'LineWidth',1.5);
% colorbar
hold on
plot(real(Z_circ),imag(Z_circ),'r','LineStyle','--','LineWidth',1.5)
scatter(real(exp(1i*pi*sin(DoA_rad))),imag(exp(1i*pi*sin(DoA_rad))),125,'r','x','LineWidth',2)
h=get(gca,'Children'); % grab all the axes handles at once
legendstr={'Null Spectrum','Unit Circle','True DOA'};
legend(h([2 1]),legendstr{[2 3]},'location','southwest')

hold off
axis equal
colormap parula
xlabel('Real','FontSize',14)
ylabel('Imaginary','FontSize',14)
zlabel('Magnitude (dB)','FontSize',14)
str = strcat( 'Null Spectrum');
title(str,'FontSize', 18)
text(-1.9,1.25, 'a)', 'FontSize', 18)
grid on
% text(-.25,.84,'1','interpreter','latex','FontSize',F_Size)
% text(-.25,-.84,'$1^*$', 'interpreter','latex','FontSize',F_Size)



%% WITH NOISE
clear; clc; 
rng_seed = 0;
rng(rng_seed)
PLOT_XENAKI = 0;

L                       = 3;                       %number of samples
M                       = 21;                       %sensors
K                       = 3;                        %number of signals
SNR                     = 0;                       %dB
SRC_RANGE               = 1;
RANDOM_SENSORS          = 0;

%% Generate necessary arrays
sin_thetas              = sort(mygen_DOAs(K,1/M)*2-1);  %sin(K random DoAs)
sin_thetas              = [-.126,.275,.67];
DoA_rad                 = asin(sin_thetas);           %angles, radians
DoA_deg                 = DoA_rad/(pi)*180;               %angles, degrees
DoA_phasors             = exp(1i*pi*sin_thetas);        %complex angles
if RANDOM_SENSORS
    r                       = [1:M-1]+1*rand(1,M-1);            %generate random sensor spacing
    r                       = [0,r]';                           %first sensor position is zero
else
    r = [0:M-1]';
end

%% Generate Signals
[ Y, Y2, noise, A, A2,X]= gen_signals_SNR( L,M,K,SNR,DoA_phasors,r,SRC_RANGE );
Y2                      = Y2+noise;
Ryy                     = Y2*Y2';
[U,~,~]                 = svd(Ryy);
En                      = U(:,K+1:end);
G                       = En*En';
% fun                     = @(z) (null_spec_polynomial( G, r, exp(1i*pi*sin(angle(z))) )); 
fun                     = @(z) (null_spec_polynomial( G, r, z )); 

%% Evaluate over z plane
xy_min                  = -1.5;
xy_max                  = abs(xy_min);
samples                 = 1000;                              %number of samples
f                       = linspace(-.5,.5,samples);         %coarse grid of frequencies between +/-.5
x                       = xy_min:.01:xy_max;
y                       = xy_min:.01:xy_max;
[X,Y]                   = meshgrid(x,y);
Z                       = X+(1i*Y);
Z_circ                  = exp(1i*2*pi.*f);
S_circ                  = zeros(1,length(Z_circ));
for i = 1:length(Z_circ)
    S_circ(i)                  = abs(fun(Z_circ(i)));
end
S                       = zeros(size(Z)); 
for i = 1:length(x)
    for j = 1:length(y)
        S(i,j)          = 10*log10(abs(fun(Z(i,j))));
    end
end
lsamp                   = 110;
levels                  = linspace(-100,10*log10(4*M),lsamp);
levels                  = 10.^(levels./10);
figure(1)
subplot(1,2,2)
[C,h] =contour(X,Y,S,levels,'LineWidth',1.5);
colorbar
hold on
plot(real(Z_circ),imag(Z_circ),'r','LineStyle','--','LineWidth',1.5)
hold on
scatter(real(exp(1i*pi*sin(DoA_rad))),imag(exp(1i*pi*sin(DoA_rad))),125,'r','x','LineWidth',2)
% legend('Null Spectrum','Unit circle','True DoA','Location','southwest')
axis equal
colormap parula
xlabel('Real','FontSize',14)
% ylabel('Imaginary','FontSize',14)
zlabel('Magnitude (dB)','FontSize',14)
str = strcat('SNR = ',{' '},num2str(SNR),' dB');
title(str,'FontSize', 18)
F_Size  = 16;
text(-1.9,1.25, 'b)', 'FontSize', 18)
grid on
% text(-.25,.84,'1','interpreter','latex','FontSize',F_Size)
% text(-.25,-.84,'$1^*$', 'interpreter','latex','FontSize',F_Size)



