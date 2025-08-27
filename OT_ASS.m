clc
clear 
close all

%%%% M.BILAL RAZA %%%% 
%%%% 22-MS-EE-07  %%%%

% Define signal parameters

K=64;   % Sampling Frequency
C=1;    % No. of Cycles
F=50;   % Fundamental Frequency
Ts=1/(K*F); % Sampling Time
k=0:C*(K-1); % Samples for Graph Plotting
% Signal To Noise Ratio in db
%SNR=0;
SNR=20;
ActualSig=zeros(1,length(K));  % Initialized the Actual Signal with zeros
% Generation of random integer between -5 and 5
X0= randi([-5 5],1,10);
omega=1;

% The signal from the Data Set
%      |  |  |      |      |
%      |HO|FR| AMP  | Phase|

Data = [1  50  1.5   80
        3  150 0.5   60
        5  250 0.2   45
        7  350 0.15  36
        11 550 0.1   30];
    
%%% Declaring the Actual Signal

   for i=1:size(Data,1)
      ActualSig = ActualSig + ( Data(i,3)*(sin(2*pi*k*Ts*Data(i,2)+Data(i,4)*pi/180)) );
   end
  
%%% Declaring the Noisy Signal   
 NoisySig=awgn(ActualSig,SNR,'measured')
 
%%% Declaring the Estimated Signal

   for i=1:size(Data,1)
      
       H(2*i-1,:) = sin(Data(i,2)*2*pi*k*Ts);
       H(2*i  ,:) = cos(Data(i,2)*2*pi*k*Ts);
       
   end
   
%%% Function to Perform Uncounstrainted Optimization
  [X,fval]=fminunc(@(X) fitnessFunction(X,H,NoisySig),X0)

%%% Estimation of Signal  
   EstimateSig = X*H;

%%% Plotting of the Signal
   plot(k,NoisySig,'r-');
   hold on
   plot(k,EstimateSig,'b.');
   hold on
   title('Signal Harmonic Graph (C-2) By 22MS-EE-07')
   plot(k,ActualSig,'g*');
   
%%% Labelling the Graph
   xlabel('Time(s)')
   ylabel('Amplitudes')
   legend('ActualSig','EstimateSig','NoisySig')

%%% Amplitude and Phase Calculation
Amplitude= zeros(1, 5);
Angle = zeros(1, 5);

   for n=1:5
    Amplitude(n) = sqrt(X(2*n - 1)^2 + X(2*n)^2);
    Angle(n) = (180/pi)*atan2(X(2*n),X(2*n - 1));
   end
   
Amplitude
Angle

%%% Mean Square Error
 squared_diff = (ActualSig - EstimateSig).^2;
    mse = 1/K * sum(squared_diff)
%%% Perforamnce Index
PER = sum((ActualSig - EstimateSig).^2) / sum(ActualSig.^2)
%%% Weight Least Square
WLS = omega*sum((NoisySig - EstimateSig).^2)

%%% Function for Overall Estimation of Harmonic Problem
function f = fitnessFunction(X,H,NoisySig)
    EstimateSig = X*H;
    f = sum((NoisySig - EstimateSig).^2);
end