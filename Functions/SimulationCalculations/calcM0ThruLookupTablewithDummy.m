%% Lookup table approach to MTsat with arbitrary readouts. 
function M0 = calcM0ThruLookupTablewithDummy(input_img, b1, T1, mask, echoSpacing, numReadExcitation, TR, flipAngle,dummyEchoes)

% Input T1map, noMT-weighted image, B1 map
% important change that you have numExcitation and DummyEcho
% numExcitation == full turbofactor (with number of read lines = numExcitation - DummyEcho)


% Make sure the T1, echoSpacing and Time delay (TD) are all in the same
% units ( milliseconds).
% Flip is in degrees


% inputMTw_img = GRE_sig(:,x,y)';
% b1 = one_mat; %b1;
% T1 = one_mat .*1/R1app(x,y) * 1000;
% mask =  one_mat; 
% M0 = one_mat.*Aapp(x,y);
% echoSpacing = Params.echoSpacing * 1000;
% numExcitation = Params.numExcitation;
% TR =  Params.TR * 1000;
% flipAngle = Params.flipAngle;
% dummyEchoes = Params.DummyEcho;

%% Handling for missing values
if isempty(mask)
    mask = ones(size(input_img));
end

if isempty(b1)
    b1 = ones(size(input_img));
end

%% Equations
flip_a = (flipAngle.*b1) * pi / 180; % correct for B1 and convert to radians

%% Following readout you magnetization (M2) = A1 + M1 * B , derivation at the bottom of function. 
TD = TR - (echoSpacing * (numReadExcitation));

x = cos(flip_a) ;
y = exp(-echoSpacing./T1);

B1 = (x.*y).^numReadExcitation;
% A1 = 0;
% for i = 1:numExcitation
%    A1 = A1 + M0*(x*y)^(i-1);
%    A1 = A1 - M0*(x*y)^(i)/x ;
% end

%% redo based on result from Munsch et al 2021 ihMT paper:
A1 = (1-y).* ( (1- x.^numReadExcitation .* exp(-numReadExcitation.*echoSpacing./T1)) ./ (1- x.*y) );

%% You then have some time TD for T1 relaxation, M3 = A2 + B2*M2
A2 = (1 - exp(-TD ./ T1));
B2 =  exp( -TD  ./ T1    );

%% You might have some dummy echoes played out before long echo train. M4 = A3 + M3 * B3
A3 = (1-y).* ( (1- x.^dummyEchoes .* exp(-dummyEchoes.*echoSpacing./T1)) ./ (1- x.*y) );
B3 = (x.*y).^dummyEchoes;


%% Since M4 = M1 and we want to solve for M1
%M = M0*(A3 + A2*B3 + A1*B2*B3) / (1-B1*B2*B3); % + A1*B2*B3

M0 = (input_img./sin(flip_a)) .*(1-B1.*B2.*B3)./(A3 + A2.*B3 + A1.*B2.*B3); % + A1*B2*B3

M0 = M0.*mask;






