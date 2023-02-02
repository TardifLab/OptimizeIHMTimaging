function [fitData1s, rawProc] = CR_loadAndCleanData( fn1, fn2)
%     
%     fn1='E:\GitHub\Bloch_simulation_Code\fitTissParams\MTsat_vals2fit.mat';
%     fn2='E:\GitHub\Bloch_simulation_Code\fitTissParams\MTsat_vals2fit_Mar16.mat';

d1 = load(fn1);
d2 = load(fn2);

d1 = d1.exportMat;
d2 = d2.exportMat2;

d = [d1,d2];

%% Sort data by B1;
[~, sortidx_b1] = sort(  d(10,:) , 'ascend');

% Sort the whole matrix:
mat_sort = d( :, sortidx_b1);

% Smooth matrix
mat_ss = smoothdata( mat_sort,2, 'movmedian', 10);

rawProc = mat_ss;
% figure; heatscatter(rawProc(10,:)', rawProc(5,:)'); %ylim([0 0.04]); xlim([0 15])

%% Let round B1 values to the 0.05 and then group.

mat_ss(10,:) = round( mat_ss(10,:),2);

minB1 = round(min(mat_ss(10,:)), 2);
maxB1 = round(max(mat_ss(10,:)), 2);
edges = minB1:0.01:maxB1;

% go through and group
N = histcounts(mat_ss(10,:),edges); 

% lets say anything under 1000 is noise, so we will remove.
for i = 1:length(N)

    if N(i) < 3500
        b1L = edges(i);
        b1H = edges(i+1);

        mat_ss(:,(mat_ss(10,:) >= b1L & mat_ss(10,:) < b1H  )) = [];
    end
end

%% We really just want one value for the vector inner products, so lets average. 
minB1 = round(min(mat_ss(10,:)), 2);
maxB1 = round(max(mat_ss(10,:)), 2);
edges = minB1:0.01:maxB1;

% go through and group
N = histcounts(mat_ss(10,:),edges); 
fitData1 = zeros(11, length(N));

for i = 1:length(N)
 
    b1L = edges(i);
    b1H = edges(i+1);

    temp = mat_ss(:,(mat_ss(10,:) >= b1L & mat_ss(10,:) < b1H  ));
    fitData1(:,i) = mean(temp,2);
end

% for zero cases, we get NaN, remove
fitData1(:, all(isnan(fitData1),1)) = [];
fitData1s = smoothdata( fitData1,2);
fitData1s(10,:) = fitData1(10,:);

% remove anything with really low B1:
rawProc(:,rawProc(10,:) < 0.75) = [];
rawProc(:,rawProc(10,:) > 1.1) = [];

% viewIdx = 3;
% figure; heatscatter(rawProc(10,:)', rawProc(viewIdx,:)'); %ylim([0 0.04]); xlim([0 15])
% hold on
% plot(fitData1(10,:)', fitData1(viewIdx,:)', "Color",'b',"LineWidth",3); 
% plot(fitData1s(10,:)', fitData1s(viewIdx,:)', "Color",'r',"LineWidth",2); 












