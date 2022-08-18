function [outSig, M, time_vect] = BlochSim_MPRAGESequence(Params, varargin)

%% Overview of how the code works with applicable function calls:

% Sequence timing:
% TI = time between 180 inversion pulse, and middle of excitation train(TF/2)
% ET = evolution time, between 180 inversion and first excitation
% EBT = excitation block timing (numExcitation * echospacing) 
% TD = time delay after last excitation pulse (and echospacing) until next
%   inversion
% TR = ET + EBT + TD

% if centric, dummy echoes == 2, else ==0. 

%% Note we assume that the Params.numExcitation does not include dummy scans
%% at beginning of train, but will include any at the end!

%% Need to define defaults for:
% Params.InvPulseDur
% Params.InversionEfficiency

%% If you don't want to look at the impact of the bound pool, set 
%% Params.M0b = 0, then Params.Ra == Params.Raobs;


% Params.CalcVector = 1;

%% Use name-value pairs to override other variables set. Great for parfor loops!
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        Params.(varargin{i}) = varargin{i+1};
    end
end

if ~isfield(Params,'IncludeDipolar')
    Params.IncludeDipolar = 1; 
end

if ~isfield(Params,'InversionEfficiency')
    Params.InversionEfficiency = 0.96; 
end

if ~isfield(Params,'CalcVector')
    Params.CalcVector = 0;
end


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


% if centric, dummy echoes == 2, else ==0. 
if strcmp(Params.Readout, 'centric')
    Params.DummyEcho = 2;
else % assume linear encoding
    Params.DummyEcho = 0;
end
Params.numExcitation = Params.numExcitation + Params.DummyEcho;

%% Timing Variables
% excitation block timing (Turbofactor * echospacing) 
EBT = (Params.numExcitation)*( Params.echoSpacing);

% Inversion time needs to be specified by user
TI = Params.TI; 

% evolution time, between 180 inversion and first excitation
if strcmp(Params.Readout, 'centric')
    ET = TI ; 
else % assume linear encoding
    ET = TI - EBT/2; 
end

% time delay after last excitation pulse (and echospacing) until next inversion
TD = Params.TR - ET - EBT;

if TD < 0 
    error('Check timing, calculated TD < 0')
end

% Impact of Excitation Pulse of Bound pool
Params = CalcBoundSatFromExcitationPulse(Params, Params.flipAngle); % Rrfb_exc and Rrfd_exc for excitation pulses
Params = CalcBoundSatFromInversionPulse(Params); % Rrfb_exc and Rrfd_exc for inversion pulses

%% Need to determine a sufficient number of isochromats. 
% can use this function to use more if needed.
% Params.N_spin = DetermineNumberIsoChromat(Params, TD)

if Params.PerfectSpoiling % number of spins wont matter in this case
    Params.N_spin = 1;
else
    Params.N_spin = 201;
end


%% Setup Matrices
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
    
    %% Start by applying an inversion pulse

    R = RotationMatrix_withBoundPool_Inversion( pi*Params.InversionEfficiency, 0, Params);

    % Instanteous RF pulse
    M_t = pagemtimes(R,M_t); % 50 percent faster than loop

    % for ns = 1:N_spin
    %    M_t(:,ns) = R*squeeze(M_t(:,ns));
    % end

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1);
        idx = idx+1;
    end % For viewing; 

    %% Spin evolution over time ET
    M_t = XYmag_Spoil(Params, M_t, ET, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1)+TD;
        idx = idx+1;
    end % For viewing; 


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
        % Note that Params.echoSpace ~= 0.
        M_t = XYmag_Spoil( Params, M_t, Params.echoSpacing, 0, 1);
            
        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1)+ Params.echoSpacing;
            idx = idx+1;
        end % For viewing;      
        
    end % End '1: Params.numExcitation' 

    %% Spin evolution over time TD
    M_t = XYmag_Spoil(Params, M_t, TD, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1)+TD;
        idx = idx+1;
    end % For viewing; 

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



