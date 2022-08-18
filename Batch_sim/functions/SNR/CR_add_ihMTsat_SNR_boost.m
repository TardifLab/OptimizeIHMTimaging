function output = CR_add_ihMTsat_SNR_boost(simResultsB, ParameterSet, noiseLvl, T1, MP2RAGE)
% breaks up the simResultsB becuase my computer keeps crashing when
% running the whole thing. Do in 10 steps, save each step, then combine for
% the output. 

DATADIR = pwd;
DATADIR = [DATADIR,'/intermediateResults/'];
mkdir(DATADIR)

if nargin < 3 || isempty(noiseLvl)
  noiseLvl = 0.0005; 
end
if nargin < 4 || isempty(T1)
  T1 = 1.4;
end
if nargin < 5 || isempty(MP2RAGE)
    MP2RAGE.B0=3;           % in Tesla
    MP2RAGE.TR=5;           % MP2RAGE TR in seconds 
    MP2RAGE.TRFLASH=6.4e-3; % TR of the GRE readout
    MP2RAGE.TIs=[940e-3 2830e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.NZslices=176;% C.R. just need the turbofactor here.
    MP2RAGE.FlipDegrees=[4 5];% Flip angle of the two readouts in degrees
end


divNum = floor(size(simResultsB,1)/20);
output = zeros(size(simResultsB,1), size(simResultsB,2)+2);


sIdx = 1;
eIdx = divNum;
tic
for i = 1:20

    temp = CR_add_ihMTsat_SNR(simResultsB(sIdx:eIdx,:),ParameterSet(sIdx:eIdx,:), noiseLvl, T1, MP2RAGE);

    output(sIdx:eIdx,:) = temp; % fill output matrix

    save(strcat(DATADIR,'simResults_3T_batch_boost_withCalc_intermed',num2str(i),'.mat'),'output' );
    disp( [ num2str(i/20*100),'% done boosted sims'])
    toc

    if i < 20
        sIdx = eIdx + 1;
        eIdx = eIdx+ divNum;
    
        if i == 19
            eIdx = size(simResultsB,1); % for last iteration run the rest
        end
    end

end

