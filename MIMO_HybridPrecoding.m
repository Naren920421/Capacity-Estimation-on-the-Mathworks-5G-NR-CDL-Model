function Capacity = MIMO_HybridPrecoding(H,At,Ar,BW,TX_power_mmWave,Noise_power,Ns)
%%% This is the function to calculate the channel capacity for millimeter wave signals.
%%% The hybrid precoding algorithms are proposed in [1].
% 
% Input
%           H                   : narrowband channel matrix, Nr x Nt
%           At/Ar               : the collection matrix of steering vectors of BS/UE, Nt/Nr x MN
%           BW                  : Hz
%           TX_power_mmWave     : dBm
%           Noise_power         : dBm/Hz
%           Ns                  : number of data streams, 1
% Ouput
%           Capacity            : capacity in bits/s
%
% [1] O. E. Ayach et al.,"Spatially sparse precoding in millimeter wave MIMO systems," Mar. 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of RF chains = Number of data streams
NumRF = Ns;

Nr = size(Ar,1);

TX_power_mmWave = 10^(TX_power_mmWave/10)*1e-3; % watt
Noise_power = 10^(0.1*Noise_power)*1e-3; % watt/Hz

%%%  Optimal precoder

% Singular value decomposition of the channel
[~,~,V] = svd(H);

% Optimal precoder
F_Opt = V(:,1:Ns); % Nt x Ns    
    
%%% Hybrid precoding at the BS (Algorithm 1)       

F_RF = [];
F_Res = F_Opt;
for r = 1:NumRF             
    Psi = At'*F_Res; % MN x Ns
    [~,k] = max(diag(Psi*Psi')); % MN x MN
    F_RF = [F_RF At(:,k)]; 
    F_BB = (F_RF'*F_RF)\(F_RF'*F_Opt); 
    F_Res = (F_Opt-F_RF*F_BB)/norm(F_Opt-F_RF*F_BB,'fro'); 
end       
    F_BB = sqrt(Ns)*(F_BB/norm(F_RF*F_BB,'fro')); 
        
        
%%% Hybrid combining at the UE (Algorithm 2)

% Covariance matrix of the receive symbols at the receive antennas
CovRx = (TX_power_mmWave/Ns)*H*F_RF*F_BB*F_BB'*F_RF'*H'+Noise_power*eye(Nr); % Nr x Nr

% Optimal unconstrained MMSE combiner without the hardware limitation
W_MMSE = ((1/sqrt(TX_power_mmWave))*(F_BB'*F_RF'*H'*H*F_RF*F_BB+((Noise_power*Ns)/TX_power_mmWave)*eye(Ns))\(F_BB'*F_RF'*H'))';

W_RF = [];
W_Res = W_MMSE;        
for r = 1:NumRF
    Psi = Ar'*CovRx*W_Res; 
    [~,k] = max(diag(Psi*Psi'));
    W_RF = [W_RF Ar(:,k)];
    W_BB = (W_RF'*CovRx*W_RF)\(W_RF'*CovRx*W_MMSE);
    W_Res = (W_MMSE-W_RF*W_BB)/norm(W_MMSE-W_RF*W_BB,'fro');
end

%%% Capacity calculation

% Noise pwoer across the whole bandwidth
Rn = BW*Noise_power*W_BB'*W_RF'*W_RF*W_BB;  

% Capacity in bits/s
Capacity = BW*abs(log2(det(eye(Ns)+(TX_power_mmWave/Ns)*(Rn\(W_BB'*W_RF'*H*F_RF*F_BB*F_BB'*F_RF'*H'*W_RF*W_BB)))));


    
        





