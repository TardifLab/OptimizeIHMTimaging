function [xVect1, yVect1, zEff1, zEff2, zAbs1, zAbs2] = ...
    ihMTsim_SNRvectorGen4Figures( ParameterSet , idx, zidx, simResults )

% this function takes your inputs and extracts the max values for plotting
% a grid

% Inputs
% ParameterSet = ParameterSet for sequences
% idx = the column number of interest to sort the ParameterSet (pick 2)
% zidx = 4 number vector that contains the column numbers/index for ...
%             [ihMT proxy, ihMT proxy efficiency, ihMTsat, ihMTsat
%             efficiency]  (example 11:14)


% Outputs
% xVect is the output vector sorted for x dim
% yVect is the output vector sorted for y dim
% zEff1 is the contrast efficiency sorted to match xVect and yVect
% zAbs1 is the contrast efficiency sorted to match xVect and yVect 
% zEff2 is the ihMTsat efficiency sorted to match xVect and yVect
% zAbs2 is the ihMTsat efficiency sorted to match xVect and yVect


x = ParameterSet(:,idx(1));
y = ParameterSet(:,idx(2));

comb_mat = [ x, y, simResults];

xu = unique(x);
yu = unique(y);

zidx = zidx +2; % account for concatenating x, and y

% find the max value (column 1), matrix idx (column 2)
ihMTRsnr_maxAbs = zeros( length(xu), length(yu) );
ihMTRsnr_maxEff = zeros( length(xu), length(yu));
ihMTsatsnr_maxAbs = zeros( length(xu), length(yu) );
ihMTsatsnr_maxEff = zeros( length(xu), length(yu));
xVect1 = zeros( length(xu), length(yu) );
yVect1 = zeros( length(xu), length(yu));

for i = 1:length(xu) % stack x in the 3rd dimension
   
    temp1 = comb_mat( x== xu(i),: );
    
    for j = 1:length(yu) % for each simulated number of excitations per offset Resonance
        
        temp2 = temp1(temp1(:,2)== yu(j),: );
        
        if ~isempty(temp2) % incorparating time restrictions makes this necessary
            xVect1( i, j ) = xu(i);
            yVect1( i, j ) = yu(j);
            
            ihMTRsnr_maxAbs( i, j ) = max(temp2( :, zidx(1) )); % max ihMT for selected excitation number       
            ihMTRsnr_maxEff( i, j )  = max(temp2( :, zidx(2) )); % max ihMT for selected excitation number  
            ihMTsatsnr_maxAbs( i, j ) = max(temp2( :, zidx(3) )); % max ihMT for selected excitation number       
            ihMTsatsnr_maxEff( i, j )  = max(temp2( :, zidx(4) )); % max ihMT for selected excitation number  
        end
    end    
end


%% Convert to vector and remove 0's

xVect1 = xVect1(:);
yVect1 = yVect1(:);
zEff1 = ihMTRsnr_maxEff(:);
zEff2 = ihMTsatsnr_maxEff(:);
zAbs1 = ihMTRsnr_maxAbs(:);
zAbs2 = ihMTsatsnr_maxAbs(:);

% remove zeros 
xVect1( zEff1 == 0) = [];
yVect1( zEff1 == 0) = [];
zAbs2( zEff1 == 0) = [];
zAbs1( zEff1 == 0) = [];
zEff2( zEff1 == 0) = [];
zEff1( zEff1 == 0) = [];













