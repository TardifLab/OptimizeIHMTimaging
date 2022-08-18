%% Lookup table approach to MTsat with arbitrary readouts. 
function MTsat = calcMTsatThruLookupTablewithDummyV3(inputMTw_img, b1, T1, mask, M0, echoSpacing, numExcitation, TR, flipAngle,dummyEchoes)

% New for V3 -> getting consistency with numExcitation with MAMT_model_2007_7
% important change that you have numExcitation and DummyEcho
% numExcitation == full turbofactor (with number of read lines = numExcitation - DummyEcho)


% Make sure the T1, echospacing and Time delay (TD) are all in the same
% units ( milliseconds).
% Flip is in degrees


% inputMTw_img = GRE_sig(:,x,y)';
% b1 = one_mat; %b1;
% T1 = one_mat .*1/R1app(x,y) * 1000;
% mask =  one_mat; 
% M0 = one_mat.*Aapp(x,y);
% echoSpacing = Params.echospacing * 1000;
% numExcitation = Params.numExcitation;
% TR =  Params.TR * 1000;
% flipAngle = Params.flipAngle;
% dummyEchoes = Params.DummyEcho;




%% Handling for missing values
if isempty(mask)
    mask = ones(size(inputMTw_img));
end

if isempty(b1)
    b1 = ones(size(inputMTw_img));
end

%% First divide by the MTw image by M0 to get relative signal 
MTw_sig = inputMTw_img ./ M0;

TD = TR - (echoSpacing * (numExcitation));

%% Generate a lookup table based on B1, T1 and MTsaturation

%setup vectors
B1_vector = 0.005:0.025:1.9; % get B1 contour map style artifact in sat maps with higher increment
T1_vector = (0.5:0.015:5) *1000; % Increment from 0.025 to 0.015 produced tissue dependent changes on the order of 5e-6
MT_vector = 0:0.005:1; % Increment from 0.005 to 0.001 produced tissue dependent changes on the order of 2e-6
MTsig_vector = 0:0.0005:0.25; % Increase end from 0.25 to 0.35 and (separately) increment from 0.0005 to 0.0001 got changes on order of 3e-6

SignalMatrix = zeros(length(MT_vector),1);
MTsatMatrix = zeros(length(B1_vector), length(T1_vector), length(MTsig_vector));


% calculate the lookup table
for i = 1:length(B1_vector)
    
    for j = 1:length(T1_vector)
        
        for k = 1:length(MT_vector)
            
            SignalMatrix(k) = MTrage_sig_eqn_withDummy(echoSpacing, flipAngle, T1_vector(j), TD, numExcitation - dummyEchoes, 1, MT_vector(k), B1_vector(i), 1, dummyEchoes);  

        end
        
        % We have the signal values for the metrics. We need the total
        % matrix to have the z direction be signal with the matrix values
        % being the MTdrop value. 
        MTsatMatrix(i,j,:) = interp1(SignalMatrix , MT_vector, MTsig_vector, 'pchip',0);     
    end
end

%% Now fit the image using gridded interpolant
% matrix values (MTsat) are defined by vectors: B1, T1 and MTw signal
[b, t, m] = ndgrid(B1_vector, T1_vector, MTsig_vector);
F = griddedInterpolant(b ,t, m, MTsatMatrix);

%% Turn the images into vectors then fit
q = find( (mask(:)>0));
b1_v = b1(q);
t1_v = T1(q);
mt_v = MTw_sig(q);


mtsat = F(b1_v, t1_v, mt_v);

MTsat = zeros( size(T1));
MTsat(q) = mtsat;




