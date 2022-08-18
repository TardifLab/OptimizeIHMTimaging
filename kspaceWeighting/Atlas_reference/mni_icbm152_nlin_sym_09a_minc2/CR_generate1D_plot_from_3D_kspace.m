function CR_generate1D_plot_from_3D_kspace(inputKspace, dimension)

% inputKspace = fftshift(log( abs(wm_fft)) );

if ~exist('dimension', 'var')
    dimension = 1;
end

if ~isreal(inputKspace)
    inputKspace = abs(inputKspace);
end

% need to do something special for 3D
[x, y, z] = ind2sub(size(inputKspace),find(inputKspace == max(inputKspace, [], 'all')));


if dimension == 1
    
    vp = squeeze(inputKspace(:,y,z));
    
elseif dimension == 2
    vp = squeeze(inputKspace(x,:,z));
        
else
    vp = squeeze(inputKspace(x,y,:));
        
end

figure;
plot(vp, 'LineWidth', 2)










