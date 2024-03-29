function [inputMag, satFlipAngle] = CR_batch_simSequenceFunction(Params, varargin)
%% A couple of edge cases:
% Not enough time: b1 set to 0, inputMat set to 0;
% Excitation pulses use all the SAR: b1 set to -1, inputMat set to 0;


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

B1peak_limit = 28e-6; % in Tesla

if (strcmp( Params.SatPulseShape, 'gaussian'))
    Params.PulseOpt.bw = 2./Params.pulseDur;
elseif( strcmp(Params.SatPulseShape ,'gausshann'))
    Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
end

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
     satFlipAngle = 0; % no MT pulse
     
      if necessTime > Params.TR
            inputMag = 0;           
            return;
      else
          Params.boosted = 0;
          [inputMag, ~, ~] = BlochSimFlashSequence_v2(Params,'satFlipAngle',0,'MTC',0); % reference signal simulation
      end


  %% If boosted protocol   
 elseif Params.boosted

    % time of sat pulses
    satTime = (Params.numSatPulse*(Params.pulseDur+Params.pulseGapDur));

    if satTime > Params.TR_MT % if sat scheme doesn't fit, go to next value...
        inputMag = 0;
        satFlipAngle = 0;
        return;
    else

        % If sat time fits, check whole TR
        % (remove the last TD_MT)
        Params.TD_MT = Params.TR_MT - Params.numSatPulse* (Params.pulseDur + Params.pulseGapDur) ;
        necessTime = (Params.TR_MT*Params.satTrainPerBoost)-Params.TD_MT + Params.numExcitation*(Params.echoSpacing); 

        if necessTime > Params.TR
            inputMag = 0;
            satFlipAngle = 0;
            return;
        else

            % Calculate B1rms
            Params = CR_SAR_scale_PulseHeight(Params, B1peak_limit);

            satFlipAngle = Params.satFlipAngle;

            if (satFlipAngle < 0) 
                inputMag = 0;
            else
                [inputMag, ~, ~] = BlochSimFlashSequence_v2(Params); % MT-weighted signal simulation    
            end
                
            return;
        end
     end % End else for necessary Time


%% For the non-boosted protocol                            
else

    % check timing, return if it doesnt fit
    necessTime = Params.numSatPulse*(Params.pulseDur+Params.pulseGapDur)+ Params.numExcitation*(Params.echoSpacing); 

    if necessTime > Params.TR
        inputMag = 0;
        satFlipAngle = 0;
        return;
    else

        % Calculate B1rms
        Params = CR_SAR_scale_PulseHeight(Params, B1peak_limit);

        satFlipAngle = Params.satFlipAngle;

        if (satFlipAngle < 0) 
                inputMag = 0;
            else
                [inputMag, ~, ~] = BlochSimFlashSequence_v2(Params); % MT-weighted signal simulation    
        end

        return;

    end % end if necessaryTime
end % end if Params.boosted
                        
