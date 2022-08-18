function [inputMag, B1_val] = CR_batch_simSequenceFunction(Params, varargin)

% 
% Params.delta = ParameterSet(qi,1);
% Params.flipAngle= ParameterSet(qi,2);
% Params.TR= ParameterSet(qi,3);
% Params.numSatPulse= ParameterSet(qi,4);
% Params.pulseDur= ParameterSet(qi,5);
% Params.numExcitation= ParameterSet(qi,6);
% Params.satTrainPerBoost= ParameterSet(qi,7);
% Params.TR_MT= ParameterSet(qi,8);
% Params.freqPattern='single';


%% Use name-value pairs to override other variables set. Great for parfor loops!
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        Params.(varargin{i}) = varargin{i+1};
    end
end

B1rms_limit = 14e-6; % in Tesla

Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.

%% Dummy echoes the way the TFL sequence adds them
if Params.numExcitation == 1
    Params.DummyEcho = 0;
elseif Params.numExcitation < 4
    Params.DummyEcho = 1;
else 
    Params.DummyEcho = 2;
end

Params.numExcitation = Params.numExcitation + Params.DummyEcho;

if ~isfield(Params,'ReferenceScan') % Assume Sat Pulses are played out, unless specified otherwise.
    Params.ReferenceScan = 0; 
end



 %% If reference scan for MTR calc
 if Params.ReferenceScan % used for calculating MTR.
     necessTime = Params.numExcitation*(Params.echoSpacing); 
     B1_val = 0; % no MT pulse
     
      if necessTime > Params.TR
            inputMag = 0;           
            return;
      else
          Params.boosted = 0;
          [inputMag, ~, ~] = BlochSimFlashSequence_v2(Params,'b1',0,'MTC',0); % reference signal simulation
      end


  %% If boosted protocol   
 elseif Params.boosted

    % time of sat pulses
    satTime = (Params.numSatPulse*(Params.pulseDur+Params.pulseGapDur));

    if satTime > Params.TR_MT % if sat scheme doesn't fit, go to next value...
        inputMag = 0;
        B1_val = 0;
        return;
    else

        % If sat time fits, check whole TR
        % (remove the last TD_MT)
        Params.TD_MT = Params.TR_MT - Params.numSatPulse* (Params.pulseDur + Params.pulseGapDur) ;
        necessTime = (Params.TR_MT*Params.satTrainPerBoost)-Params.TD_MT + Params.numExcitation*(Params.echoSpacing); 

        if necessTime > Params.TR
            inputMag = 0;
            B1_val = 0;
            return;
        else

            % Calculate B1rms
            Params = CR_SAR_scale_PulseHeight(Params);

            % check limit, reset if necessary
            if Params.satRMS > B1rms_limit
                Params.satRMS = B1rms_limit;                               
            end
            Params.b1 = Params.satRMS *1000000; % convert to microTesla
            B1_val = Params.b1;

            [inputMag, ~, ~] = BlochSimFlashSequence_v2(Params); % MT-weighted signal simulation        
            return;
        end
     end % End else for necessary Time


%% For the non-boosted protocol                            
else

    % check timing, return if it doesnt fit
    necessTime = Params.numSatPulse*(Params.pulseDur+Params.pulseGapDur)+ Params.numExcitation*(Params.echoSpacing); 

    if necessTime > Params.TR
        inputMag = 0;
        B1_val = 0;
        return;
    else

        % Calculate B1rms
        Params = CR_SAR_scale_PulseHeight(Params);

        % check limit, reset if necessary
        if Params.satRMS > B1rms_limit
            Params.satRMS = B1rms_limit;                               
        end
        Params.b1 = Params.satRMS *1000000; % convert to microTesla
        B1_val = Params.b1;


        [inputMag, ~, ~] = BlochSimFlashSequence_v2(Params); % MT-weighted signal simulation
        return;

    end % end if necessaryTime
end % end if Params.boosted
                        
                        
                        
                        
                        
