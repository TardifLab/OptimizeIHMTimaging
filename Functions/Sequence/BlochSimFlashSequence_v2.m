function [outSig, M, time_vect] = BlochSimFlashSequence_v2(Params, varargin)

% V2 carries the isochromats all the way through

%% Overview of how the code works with applicable function calls:
% 1. IF MTC -> call to Bloch_McConnell_wDipolar( Params, delta, b1) to get
%    RF saturation and spin evolution matrix

% 2. IF MTC -> call to SpinEvolution_Relaxation( Params, M_in, t);  for 
%    spin evolution for time gap

% 3. IF MTC & Last sat pulse done, call to XYmag_Spoil( Params, M_in, Params.MTSpoilTime, 1 )
%    for gradient spoiling, spin diffusion and spin evolution (option to set perfect spoil). 


% 4. Call to RotationMatrix_withBoundPool(Params.flipAngle*pi/180, RFphase(j)*pi/180, Params)
%    For instanteous excitation with a rotation matrix applied to the water
%    pool, and Bloch-McConnell style saturation of bound (and dipolar) pool.
%    This includes RF spoiling. 

% 5. Call to XYmag_Spoil( Params,M_in, Params.MTSpoilTime, 1 )
%    for gradient spoiling, spin diffusion and spin evolution (option to set perfect spoil). 

% 6. Call to SpinEvolution_Relaxation( Params, M_in, t);  for 
%    spin evolution for time gap


%% Use name-value pairs to override other variables set. Great for parfor loops!
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        Params.(varargin{i}) = varargin{i+1};
    end
end

if ~isfield(Params,'IncludeDipolar')
    Params.IncludeDipolar = 1; 
end


% if ~isfield(Params,'kf')
%     Params.kf = (Params.R*Params.M0b); 
% end
% 
% if ~isfield(Params,'kr')
%     Params.kr = (Params.R*Params.M0a);
% end
Params.kf = (Params.R*Params.M0b); 
Params.kr = (Params.R*Params.M0a);

if isempty(Params.Ra) % allow you to specify either Ra or Raobs
    Params.Ra = Params.Raobs - ((Params.R * Params.M0b * (Params.R1b - Params.Raobs)) / (Params.R1b - Params.Raobs + Params.R));
    if isnan(Params.Ra)
        Params.Ra = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build sequence, then convert to loop structure.
% play sequence for 5 seconds, then fill sampling table.
stepSize = 50e-6; % 50 microseconds
RFphase = 0; % starting excitation phase
last_increment = 0;
num2avgOver = 20; % you get some variation in signal, so keep the last few and average

% Equivalent of 5 seconds of imaging to steady state, then record data. 
loops = ceil(6/Params.TR) + num2avgOver;

%% Standard Stuff
M0 = [0 0 Params.M0a, Params.M0b, 0]';
I = eye(5); % identity matrix      
B = [0 0 Params.Ra*Params.M0a, Params.R1b*Params.M0b, 0]';

if Params.echoSpacing == 0
    Params.echoSpacing = 5e-3; % ensure value not 0
end    


if Params.MTC && Params.boosted
    TD = Params.TR -  (Params.TR_MT*Params.satTrainPerBoost) + Params.TD_MT - (Params.numExcitation*Params.echoSpacing) - Params.G_time_elapse_MT; % time in seconds
elseif Params.MTC
    TD = Params.TR - (Params.numSatPulse)*( Params.pulseDur + Params.pulseGapDur) - Params.numExcitation*( Params.echoSpacing) ...
         + Params.pulseGapDur - Params.G_time_elapse_MT; % time in seconds
else % with no MTC
    if Params.numExcitation == 1
        TD = 0;
        Params.echoSpacing = Params.TR; % For GRE and/or FLASH style sequence
        Params.DummyEcho = 0;
    else
        TD = Params.TR - Params.numExcitation*( Params.echoSpacing); % time in seconds
    end
end

if TD < 0
    error('Check timing variables, TD < 0');
end


%% Temporary code to hold over until fix old code:
if Params.MTC && strcmp(Params.SatPulseShape,'hanning')
    Params.SatPulseShape = 'gausshann';
    disp('Switching SatPulseShape from hanning to gausshann...')
    disp('Best to update your code to reflect switch to qMRlab base')
end


%% Precompute MTC Pulses:

if Params.MTC
    tSat = 0 : stepSize : Params.pulseDur;

    PulseDur = ceil(Params.pulseDur/stepSize); % Break down pulse into rectangles
    alpha = Params.satFlipAngle;

    if ~isfield(Params,'PulseOpt')
        Params.PulseOpt = [];
    end
    if isempty(Params.PulseOpt)
        disp('Using default sat pulse bandwidth, otherwise set Params.PulseOpt.bw to a value in Hz')
    end
    satPulse = GetPulse(alpha, Params.delta, Params.pulseDur, Params.SatPulseShape, Params.PulseOpt);

    E_rf = zeros(5,5,PulseDur);

    % Precompute the RF matrix that is time variant in the 3rd dimension
    for k = 1:PulseDur
        E_rf(:,:,k) = Bloch_McConnell_wDipolar(Params, Params.delta, satPulse.omega(tSat(k))); 
    end
    
    if strcmp( Params.freqPattern,'dualAlternate')

        E_rf2 = zeros(5,5,PulseDur); % second matrix for negative

        % Precompute the RF matrix that is time variant in the 3rd dimension        
        for k = 1:PulseDur
            E_rf2(:,:,k) = Bloch_McConnell_wDipolar(Params, -Params.delta, satPulse.omega(tSat(k))); 
        end
    end

end


% Impact of Excitation Pulse of Bound pool
Params = CalcBoundSatFromExcitationPulse(Params, Params.flipAngle); % Rrfb_exc and Rrfd_exc for excitation pulses

%% Need to determine a sufficient number of isochromats. 
% can use this function to use more if needed.
% Params.N_spin = DetermineNumberIsoChromat(Params, TD)

if Params.PerfectSpoiling % number of spins wont matter in this case
    Params.N_spin = 1;
else
    Params.N_spin = 201;
end


%% Setup Matrices
Params.CalcVector = 0;

M = zeros(5,loops*20);
M(:,1) = M0;

time_vect = zeros( loops*20,1);
M_t = repmat(M0,1, Params.N_spin);

Sig_vec = zeros(num2avgOver, Params.numExcitation-Params.DummyEcho );
rep = 1; % to count over the number to average over

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of sequence loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 2;

for i = 1:loops
    
    %% MT block:
    if Params.MTC    
    
        %% Boosted protocol saturation  
        if Params.boosted
                   
            for j1 = 1:Params.satTrainPerBoost % repeat Sat pulse train for specified number of times for boost   
                for j = 1:Params.numSatPulse % for each MTsat pulse cycle
    
                    % For Saturation Pulse
                    for k = 1:PulseDur
                        if strcmp( Params.freqPattern,'dualAlternate') && (mod(j,2) == 0)
                            [M_t, time_vect(idx)] = calcPoolChange(E_rf2(:,:,k), B, I, stepSize, M_t, time_vect(idx-1) );
                        else % single is default;   
                            [M_t, time_vect(idx)] = calcPoolChange( E_rf(:,:,k), B, I, stepSize, M_t, time_vect(idx-1) );
                        end
   
                        if Params.CalcVector == 1 
                            M(:,idx) = mean(M_t,2); 
                            time_vect(idx) = time_vect(idx-1) +  stepSize;
                            idx = idx +1;
                        end % For viewing; 
                        
                    end    
    
                    % During the Pulse Gap, the pools relax, Rrfb = 0
                    M_t = XYmag_Spoil(Params, M_t, Params.pulseGapDur, 0, 0);
                    
                    if Params.CalcVector == 1
                        M(:,idx) = mean(M_t,2); 
                        time_vect(idx) = time_vect(idx-1) +  Params.pulseGapDur;
                        idx = idx+1;
                    end
              
                end % for Params.numSatPulse
    
                if j1 < Params.satTrainPerBoost % don't include the relaxation time after last one.
               
                    % At the end of pulse train, relax some time TD_MT                   
                    %M_t = SpinEvolution_Relaxation( Params, M_t,  Params.TD_MT);
                    M_t = XYmag_Spoil(Params, M_t, Params.TD_MT, 0, 0);
                    
                    if Params.CalcVector == 1
                        M(:,idx) = mean(M_t,2); 
                        time_vect(idx) = time_vect(idx-1) +  Params.TD_MT;
                        idx = idx+1;
                    end % For viewing; 
    
                end % end if counter less than Params.satTrainPerBoost
            end % for Params.satTrainPerBoost
            
             
        %% Non-boosted protocol saturation   
        else        
            for j = 1:Params.numSatPulse % for each MTsat pulse cycle
    
                % For Saturation Pulse
                for k = 1:PulseDur
                   if strcmp( Params.freqPattern,'dualAlternate') && (mod(j,2) == 0)
                        A_sat = E_rf2(:,:,k);
                    else % single is default;   
                        A_sat = E_rf(:,:,k);
                    end
    
                    [M_t, time_vect(idx)] = calcPoolChange(A_sat, B, I, stepSize, M_t, time_vect(idx-1) );
                    
                    if Params.CalcVector == 1 
                        M(:,idx) = mean(M_t,2); 
                        time_vect(idx) = time_vect(idx-1) +  Params.pulseDur;
                        idx = idx +1;
                    end % For viewing; 
                end    
    
                % During the Pulse Gap, the pools relax, Rrfb = 0
                M_t = XYmag_Spoil(Params, M_t, Params.pulseGapDur, 0, 0);
                
                if Params.CalcVector == 1
                    M(:,idx) = mean(M_t,2);
                    time_vect(idx) = time_vect(idx-1) +  Params.pulseGapDur;
                    idx = idx+1;
                end % For viewing; 
   
            end % End 'for l:Params.numSatPulse'
        end % End 'if boosted'
    
    end % End 'if Params.MTC'    

   
    if Params.MTC
        %% Spoil - for both Boosted and non-boosted
        % If Params.PerfectSpoiling, then set XY to 0. 
        % Otherwise, compute relaxation, and spin diffusion over MTSpoilTime
        M_t = XYmag_Spoil(Params, M_t, Params.G_time_elapse_MT, 1, 1);
        if Params.CalcVector == 1 
            M(:,idx) = mean(M_t,2);  
            time_vect(idx) = time_vect(idx-1)+ Params.G_time_elapse_MT;
            idx = idx+1;
        end % For viewing; 
    end

    %% Excitation Block
    % Keep track of 5 magnetization vectors through excitation through
    % instanteous rotation of water pool, plus 'instanteous' saturation of bound pool
    % Keep track of XY mag for RF spoiling, gradient spoiling and spin diffusion. 
    % Signal == XY magnetization immediately following application of Rotation. 
    
    % Compute the RF phase for spoiling for entire excitation train
    if Params.RFspoiling
        [RFphase, last_increment] = IncrementRFspoilPhase( RFphase(end), Params, last_increment);
    else
        RFphase = zeros(1,Params.numExcitation);
    end
    
    for j = 1: Params.numExcitation
        
        % Calculate rotation matrix for excitation-specific phase
        R = RotationMatrix_withBoundPool(Params.flipAngle*pi/180, RFphase(j)*pi/180, Params);

        % Instanteous RF pulse
        M_t = pagemtimes(R,M_t); % 50 percent faster than loop

        % If you are missing pagemtimes, swap it with the slower code below:
        % for ns = 1:N_spin
        %    M_t(:,ns) = R*squeeze(M_t(:,ns));
        % end
    
        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1);
            idx = idx+1;
        end % For viewing;  

        %% Store the magnetization of each excitation pulse after 5 seconds prep
        if i > loops-num2avgOver && (j > Params.DummyEcho) 

            Sig_vec(rep,j - Params.DummyEcho ) = TransverseMagnetizationMagnitude(M_t);
    
           if (i == loops) && (j == Params.numExcitation) % if simulation is done...             
               outSig = mean(Sig_vec,1); % output 1xTurbofactor vector
               if Params.CalcVector == 1
                   M(:,idx:end) = [];
                   time_vect(idx:end) = [];
               end

               return
           end
            
           % increase repetition index
           if j == Params.numExcitation
               rep = rep+1;
           end
        end % End 'if SS_reached'   

        
        %% Apply Gradient Spoiling
        % Note that Params.echoSpace ~= 0. For Flash sequence, put all spin
        % evolution into here.
        
        M_t = XYmag_Spoil( Params, M_t, Params.echoSpacing, 0, 1);
               
        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1)+ Params.echoSpacing;
            idx = idx+1;
        end % For viewing;      
        
    end % End '1: Params.numExcitation' 
       
    %% Spin Evolution
      
    % Calculate spin evolution with diffusion
    if TD > 0
        M_t = XYmag_Spoil(Params, M_t, TD, 0, 0);

        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1)+TD;
            idx = idx+1;
        end % For viewing; 
    end


end





% %% Debug and view
% figure;
% plot(time_vect, sqrt(sum(M(1:2,:).^2)))
% % 
% figure;
% plot(time_vect, M(3,:))
% 
% figure;
% plot(time_vect, M(4,:))

% warning('') % Clear last warning message
% [warnMsg, warnId] = lastwarn;
% if ~isempty(warnMsg)
% return
% end

% 
% figure;
% plot(M_t(1,:),'-r','LineWidth',1)
% hold on
% plot(M_t(2,:),'--b')









