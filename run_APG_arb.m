%% Run Alternating Projections based Gridless beamforming for ARBITRARY array geometry
clear; clc;

%% Declare adjustable variables
trials                  = 200;                      
L                       = 1;            %number of samples
M                       = 20;           %number of sensors
K                       = 3;            %number of signals
SNR                     = 10;           %in dB
SRC_POW_RANGE           = 5;            %max amplitude difference between sources
MAX_IT                  = 10*K;         %iterations
RANDOM_SENSORS          = 1;            %turn on arbitrary array?
NOISE                   = 1;            %turn on noise in signal

%% ALGORITHMS
APG                     = 1;            %Alternating Projections based Gridless
ADMM                    = 0;            %Alternating Directions Method of Multipliers
SBL                     = 0;            %Sparse Bayesian Learning
LASSO                   = 0;            %Least Absolute Shrinkage and Selector Operator
ROOTMUSIC               = 0;            %root-MUSIC
IRREGULAR_RM            = 0;            %Irregular root-MUSIC

%% INITIALIZATIONS
T_APGs                  = zeros(trials,1);
T_RMs                   = zeros(trials,1);
T_iRMs                  = zeros(trials,1);
T_ADMMs                 = zeros(trials,1);
RMSE_iAPGs              = zeros(trials,1);
RMSE_RMs                = zeros(trials,1);
RMSE_iRMs               = zeros(trials,1);
RMSE_ADMMs              = zeros(trials,1);
RMSE_SBLs               = zeros(trials,1);                
Erec_iAPGs              = zeros(trials,1);
Erec_ADMMs              = zeros(trials,1);

%% RUN ALL
for i = 1:trials
    rng(i) %seed = i
    
    %% Generate DoAs, sensor positions, specify source power range
    sin_thetas              = sort( (gen_DOAs(K,(1/M))-.5)*2);  %sin(K random DoAs)
    DoA_phasors             = exp(1i*pi*sin_thetas);            %complex phasors
    DoA_rad                 = asin(sin_thetas);                 %angles, radians
    DoA_deg                 = DoA_rad/(pi)*180;                 %angles, degrees
    r                       = ([1:M-1])+1*(rand(1,M-1)-.5);     %generate random sensor spacing
    r                       = [0,r]';                           %first sensor position is zero
    r                       = (M-1)*(r/max(r));                 %scale to constant aperture
    if RANDOM_SENSORS == 0
        r                   = [0:M-1]';
    end
    [ Y, Y2, noise, A, A2,X]= gen_signals_SNR( L,M,K,SNR,DoA_phasors,r,SRC_POW_RANGE );
    if NOISE
        Y                   = Y+noise;                              %add noise
        Y2                  = Y2+noise;
        Y                   = Y./norm(Y);                           %normalize
        Y2                  = Y2./norm(Y2);
    end
    %% GRIDLESS METHODS
    %% AP Gridless
    if APG
        tic                                                           %run and time
        if RANDOM_SENSORS
            [T,it_eAPG]             = AP_Gridless(Y2,r,K,MAX_IT,1e-2,0);      %solve for Tu
        else
            [T,it_eAPG]             = AP_Gridless(Y,[0:M-1]',K,MAX_IT,1e-2,0);  
        end
        time_eAPG               = toc;                          %record time
        T_APGs(i)               = time_eAPG;                    %save time        
        [t_est,~]               = IVD( r, K, T );               %decompose
        DoA_est_deg             = asin(sort(t_est*2))/pi*180;   %DoA estimates in deg
        phasor_est              = exp(1i*2*pi*(t_est));         %corresponding phasors
        A_est                   = (phasor_est.').^r;            %reconstructed array steering matrix
        X_est                   = pinv(A_est)*Y2;               %reconstructed signals
        Y_hat                   = A_est*X_est;                  %reconstructed 'clean' measurements                   
        RMSE_eAPG               = norm(Y_hat-Y2,'fro')./norm(Y2,'fro');%reconstruction error
        Erec_iAPGs(i)           = RMSE_eAPG;                    %record reconstruction error
        RMSE_deg                = sqrt( mean( (DoA_deg(:)-DoA_est_deg(:) ).^2));%RMSE deg
        RMSE_iAPGs(i)           = RMSE_deg;                     %record RMSE deg
    end
    
    %% ADMM
    if ADMM
        tic
        if RANDOM_SENSORS
            Tu_ADMM             = ADMM_NUA(Y2,r,K,MAX_IT,.01);  %bare bones ADMM retrofitted for NUA
        else
            Tu_ADMM             = ADMM_gridless(Y,.001,MAX_IT,0);%run original ADMM
        end
        T_ADMMs(i)              = toc;
        [t_est,~]               = IVD( r, K, Tu_ADMM );         %decompose
        DoA_est_deg             = asin(sort(t_est*2))/pi*180;   %DoA estimates in deg
        phasor_est              = exp(1i*2*pi*(t_est));         %corresponding phasors
        A_est                   = (phasor_est.').^r;            %reconstructed array steering matrix
        X_est                   = pinv(A_est)*Y2;               %reconstructed signals
        Y_hat                   = A_est*X_est;                  %reconstructed 'clean' measurements                   
        Erec_ADMM               = norm(Y_hat-Y2,'fro')./norm(Y2,'fro');%reconstruction error
        Erec_ADMMs(i)           = Erec_ADMM;                    %record reconstruction error
        RMSE_ADMM               = sqrt( mean( (DoA_deg(:)-DoA_est_deg(:) ).^2));%RMSE deg
        RMSE_ADMMs(i)           = RMSE_ADMM;
    end
    
    %% GRIDDED METHODS
    %% SBL
    if SBL
        theta_est               = sort(do_SBLML3( Y2, K, r, 1 ));
        RMSE_SBL                = sqrt(mean( (DoA_deg(:)-theta_est(:)).^2 ))
        RMSE_SBLs(i)            = RMSE_SBL;
    end
    
    %% LASSO
    if LASSO
        theta_est               = sort(do_LASSO( Y2, r, K, .5, noise));
        RMSE_LASSO              = sqrt(mean( (DoA_deg(:)-theta_est(:)).^2 ))
    end
    
    %% SUBSPACE METHODS
    %% rootMUSIC
    if ROOTMUSIC
        tic
        [DoA_est_deg]           = rootMUSIC(Y,K);
        time_RM                 = toc;
        T_RMs(i)                = time_RM;
        RMSE_RM                 = sqrt(mean((sort(DoA_deg(:))-DoA_est_deg).^2));
        RMSE_RMs(i)             = RMSE_RM;
    end
    %%  Irregular rootMUSIC
    if IRREGULAR_RM
        tic                                                           %run and time
        [t_est_eRM]             = irregular_rootMUSIC( Y2,r,K);            %solve for Tu
        time_eRM                = toc;                                %record time
        DoA_est_deg             = t_est_eRM;                            %DoA estimates in deg
        RMSE_eRM                = sum( sqrt((DoA_deg(:)-DoA_est_deg(:) ).^2));%RMSE
        RMSE_iRMs(i)            = RMSE_eRM;
        T_iRMs(i)               = time_eRM;
    end
    
    %% Print results
    if APG
        sprintf('Trial #: %i \n eAPG_test: \n RMSE (deg): %12.8f \n Rec error: %12.8f \n Time (s):   %12.8f',i, RMSE_deg, RMSE_eAPG, time_eAPG)
    end
    if ADMM 
        sprintf('Trial #: %i \n ADMM_test: \n RMSE (deg): %12.8f \n Rec error: %12.8f \n Time (s):   %12.8f',i, RMSE_ADMM, Erec_ADMM, T_ADMMs(i))
    end
    if ROOTMUSIC
        sprintf('Trial #: %i \n root-MUSIC: \n RMSE (deg): %12.8f \n Time (s): %12.8f ',i, RMSE_RM, time_RM)
    end
    if IRREGULAR_RM
        sprintf('Trial #: %i \n eAPG_test: \n Rec error: %12.8f \n Time (s):   %12.8f',i, RMSE_eRM,time_eRM)
    end
    
    
end

if APG
    successes               = find(RMSE_iAPGs<5);
    mean_APG                = mean(RMSE_iAPGs(successes));
    APG_txt2                = sprintf('eAPG_test: \n Average RMSE (deg): %12.8f \n Average Time (s):   %12.8f \n Outliers eAPG:  %12f', mean_APG,mean(T_APGs),nnz(RMSE_iAPGs>5))
end
if SBL
    successes               = find(RMSE_SBLs<5);
    mean_SBL                = mean(RMSE_SBLs(successes));
    SBL_txt2                = sprintf('SBL_test: \n Average RMSE (deg): %12.8f  \n Outliers eAPG:  %12f', mean_SBL,nnz(RMSE_SBLs>5))
end
if ROOTMUSIC
    txt3                    = sprintf('root-MUSIC: \n Average RMSE (deg): %12.8f \n Average Time (s):   %12.8f \n Outliers eAPG:  %12f', mean(RMSE_RMs),mean(T_RMs),nnz(RMSE_RMs>5))
end
if IRREGULAR_RM
    txt4                    = sprintf('e-root-MUSIC: \n Average RMSE (deg): %12.8f \n Average Time (s):   %12.8f', mean(RMSE_iRMs),mean(T_iRMs))
end
if trials > 1
    figure()
    clf
    if APG
        hold on
        plot(RMSE_iAPGs)
        hold off
    end
    if ADMM
        hold on
        plot(RMSE_ADMMs)
        hold off
    end
    if SBL
        hold on
        plot(RMSE_SBLs)
        hold off
    end
    if ROOTMUSIC
        hold on
        plot(RMSE_RMs)
        hold off
    end
    if IRREGULAR_RM
        hold on
        plot(RMSE_iRMs)
        hold off
    end
    title('RMSE vs. Trials')
    xlabel('Trial number')
    ylabel('RMSE (deg)')
    grid on
end
