%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launcher Script
% 
% This simulation is to 
%                   1. model MIMO channels by using the Mathworks 5G NR CDL channel model in the 5G Toolbox;
%                   2. estimate the average achievable channel capacity, given different link lengths.
% 
% Simulation assumptions
%                   1. The path loss effect is included by default; the shadow fading and oxygen absorption effects can be executed as optional features [1].
%                   2. The bandwidth scales linearly with the carrier frequency.
%                   3. For millimeter wave frequencies, the hybrid precoding approach is adopted at both ends of the link [2].
% 
% [1] 3GPP, "Study on channel model for frequencies from 0.5 to 100 GHz", 3rd Generation Partnership Project (3GPP), TR 38.901, Dec.2017.
% [2] O.E.Ayach et al., "Spatially sparse precoding in millimeter wave MIMO systems", Mar.2013.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Simulation Parameters

clear variables; % Clear the workspace.

%%% Parameters for channel modelling

Fc = 30e9; % Carrier frequency in Hz
WaveLength = physconst('LightSpeed')/Fc;

BS_height = 25; % BS antenna height (macro-cell scenario)
UE_height = 1.5; % UE antenna height (outdoor UEs)

Dis2D = 50:50:300; % Horizontal BS-UE distance in meters
Dis3D = sqrt((BS_height-UE_height)^2+Dis2D.^2); % Actual BS-UE distance

%%% Parameters for capacity calculation

BW = 0.005*Fc; % Bandwith in Hz
TX_power_mmWave = 35; % Transmit power in dBm, for mmWave frequencies
TX_power_Microwave = 49; % Transmit power in dBm, for microwave frequencies
Noise_power = -174; % dBm/Hz
Ns = 1; % Number of data streams
ITER = 1000; % Number of random channel realizations

%% Additional Channel Modelling Features

% Shadow fading
ShadowFadingFlag = 0; % 0 - not included, 1 - included

% Oxygen absorption
OxygenAbsorptionFlag = 0; % 0 - not included, 1 - included

%% CDL-A Channel Model

CDL_A = nrCDLChannel;

CDL_A.DelayProfile = 'CDL-A'; 
CDL_A.DelaySpread = 200*1e-9; % 200 ns
CDL_A.CarrierFrequency = Fc;
CDL_A.ChannelFiltering = false; % For extracting channel coefficients

CDL_A.TransmitAntennaArray.Size = [4 4 1 1 1]; % [M N P Mg Ng]
CDL_A.TransmitAntennaArray.PolarizationAngles = [0,0]; % Default: [45 -45], applies when P = 2
CDL_A.ReceiveAntennaArray.Size = [2 2 1 1 1]; % [M N P Mg Ng]
CDL_A.ReceiveAntennaArray.PolarizationAngles = [0,0]; % Default: [0 90], applies when P = 2


%% Channel Parameters Extraction

cdlinfo = info(CDL_A);

%%% Antenna configurations

% BS
Nt = cdlinfo.NumTransmitAntennas;
SizeBS = CDL_A.TransmitAntennaArray.Size; % Struct
ElementSpacingBS = CDL_A.TransmitAntennaArray.ElementSpacing; % [0.5 0.5 1 1]
% Obtain antenna position vectors using Phased Array System Toolbox
ArrayBS = phased.URA('Size',SizeBS(1:2),'ElementSpacing',ElementSpacingBS(1:2)*WaveLength);
ElePosBS = getElementPosition(ArrayBS); % 3 x (Nt/P), for one polarization only
ElePosBS = repelem(ElePosBS,1,SizeBS(3)); % 3 x Nt

% UE 
Nr = cdlinfo.NumReceiveAntennas;
SizeUE = CDL_A.ReceiveAntennaArray.Size; % Struct
ElementSpacingUE = CDL_A.ReceiveAntennaArray.ElementSpacing; % [0.5 0.5 0.5 0.5]
% Obtain antenna position vectors using Phased Array System Toolbox
ArrayUE = phased.URA('Size',SizeUE(1:2),'ElementSpacing',ElementSpacingUE(1:2)*WaveLength);
ElePosUE = getElementPosition(ArrayUE); % 3 x (Nr/P), for one polarization only
ElePosUE = repelem(ElePosUE,1,SizeUE(3)); % 3 x Nr

%%% Clustered angles (AOD,AOA,ZOD,ZOA)

NumCluster = length(cdlinfo.PathDelays); % N

% Cluster angles
ClusterAOD = cdlinfo.AnglesAoD; % 1 x N
ClusterAOA = cdlinfo.AnglesAoA;
ClusterZOD = cdlinfo.AnglesZoD;
ClusterZOA = cdlinfo.AnglesZoA;

% Subpath angles
RayOffset = [0.0447 0.1413 0.2492 0.3715 0.5129 0.6797 0.8844 1.1481 1.5195 2.1551];
RayOffset = [RayOffset;-RayOffset];
RayOffset = RayOffset(:); 

NumSubPaths = length(RayOffset); % M

AngleSpreads = CDL_A.AngleSpreads;

C_ASD = AngleSpreads(1);
C_ASA = AngleSpreads(2);
C_ZSD = AngleSpreads(3);
C_ZSA = AngleSpreads(4);

% M x N
RayAOD = repmat(ClusterAOD,NumSubPaths,1)+C_ASD*repmat(RayOffset,1,NumCluster); 
RayAOA = repmat(ClusterAOA,NumSubPaths,1)+C_ASA*repmat(RayOffset,1,NumCluster); 
RayZOD = repmat(ClusterZOD,NumSubPaths,1)+C_ZSD*repmat(RayOffset,1,NumCluster); 
RayZOA = repmat(ClusterZOA,NumSubPaths,1)+C_ZSA*repmat(RayOffset,1,NumCluster); 

% Randomly couple the departure and arrival angles.
RayAOD = Coupling(RayAOD,NumSubPaths,NumCluster); 
RayAOA = Coupling(RayAOA,NumSubPaths,NumCluster);
RayZOD = Coupling(RayZOD,NumSubPaths,NumCluster);
RayZOA = Coupling(RayZOA,NumSubPaths,NumCluster);

%%% Steering vectors of depature/arrival angles

[Ar,At] = getSteeringVector(RayAOD(:),RayAOA(:),RayZOD(:),RayZOA(:),ElePosBS,ElePosUE,WaveLength);

%% Capacity Estimation

% NLOS path loss model (optional model)
[PathLoss,ShaSTD] = getPathLossNLOS(Fc*1e-9,Dis3D); % dB, vector
if ShadowFadingFlag == 1
    PathLoss = PathLoss+normrnd(0,ShaSTD);
end

Capacity = zeros(length(Dis2D),ITER);

for d = 1:length(Dis2D)
    
    for i = 1:ITER
    
        % Extract channel coefficients.
        [PathGains,~] = CDL_A();
        PathGains = squeeze(mean(PathGains,1)); % Averaging over channel snapshots
        ChannelMatrix = permute(PathGains,[3,2,1]); % Nr x Nt x N
            
        % Apply path loss.
        ChannelMatrixPL = 1/sqrt(10^(PathLoss(d)/10))*ChannelMatrix;
        
        % Oxygen absorption feature
        if OxygenAbsorptionFlag == 1
            
            OxygenLoss = Oxygen_Absorption(Fc); % dB/km
            OxygenLoss = OxygenLoss*(Dis3D(d)/1000); % dB/m
           
            ChannelMatrixFinal = 1/sqrt(10^(OxygenLoss/10))*ChannelMatrixPL;
        else
            
            ChannelMatrixFinal = ChannelMatrixPL;
            
        end
                
        % Capacity calculation (assuming narrowband channels)
        NarrowbandChannels = sum(ChannelMatrixFinal,3);
        
        if floor(Fc*1e-9) < 30 % Microwave frequencies
            
            % MIMO multiplexing
            Capacity(d,i) = MIMO_Multiplexing(NarrowbandChannels,Ns,Noise_power,BW,TX_power_Microwave); % bits/s
            
        else % Millimeter wave frequencies
            
            % Hybrid precoding
            Capacity(d,i) = MIMO_HybridPrecoding(NarrowbandChannels,At,Ar,BW,TX_power_mmWave,Noise_power,Ns); % bits/s
            
        end

    end
    
end

AveCapacity = mean(Capacity,2)*1e-6; % Mbps

%% Capacity Plot

figure();title('Average Channel Capacity on CDL-A Model');
plot(Dis2D,AveCapacity,'-o','LineWidth',2.0,'MarkerSize',8.0);
xlabel('BS-UE Horizontal Distance (meter)');
ylabel('Average Capacity (Mbps)');
grid on;



    
            
        

        




