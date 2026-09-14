lf = lf_generate_fromnyhead('montage', 'S64');


  %lf = lf_generate_fromnyhead();

 %lf = lf_generate_fromnyhead('labels', {'Nz', 'LPA', 'RPA', 'PPO9h', 'PO7', 'PO9', 'POO9h', ...
  %          'O1', 'I1', 'OI1h', 'Oz', 'Iz', 'OI2h', 'O2', 'I2', 'POO10h', ...
   %        'PO8', 'PO10', 'PPO10h', 'PPO7h', 'PO5h', 'PPO3h', 'POO5', ...
    %        'PO3h', 'PPO1h', 'POO3', 'Pz', 'POz', 'PPO2h', 'POO4',...
     %       'PO4h', 'PPO4h', 'POO6', 'PO6h', 'PPO8h', 'FFT9h', 'FC5', ...
      %      'FT7', 'FT9', 'FCC5h', 'FTT7h', 'C5', 'T7', 'CCP5h', ...
       %     'TTP7h', 'CP5', 'TP7', 'TPP7h', 'P5', 'P7', 'P9', 'P6',...
        %    'P8', 'P10', 'TPP8h', 'CP6', 'TP8', 'CCP6h', 'TTP8h',...
         %   'C6', 'T8', 'FCC6h', 'FTT8h', 'FC6', 'FT8', 'FT10', ...
          %  'FFT10h', 'FC1', 'FC3', 'FCC1h', 'FCC3h', 'C1', 'C3', ...
           % 'CCP1h', 'CCP3h', 'CPz', 'CP1', 'CP3', 'CPP1h', 'CPP3h',...
            %'CPP5h', 'P1', 'P3', 'P2', 'P4', 'CPP2h', 'CPP4h', 'CPP6h',...
%            'CP2', 'CP4', 'CCP2h', 'CCP4h', 'Cz', 'C2', 'C4', 'FCC2h', ...
 %           'FCC4h', 'FC2', 'FC4', 'FFT7h', 'FFC5h', 'F5', 'FFC3h', ...
  %          'F3', 'FFC1h', 'F1', 'FCz', 'Fz', 'FFC2h', 'F2', 'FFC4h',...
   %         'F4', 'FFC6h', 'F6', 'FFT8h', 'F7', 'AFF7h', 'AF5h', 'AF7', ...
    %        'AF3h', 'Fp1', 'AFF1h', 'AFF2h', 'AFz', 'Fpz', 'AF4h', 'Fp2',...
     %       'AF6h', 'AF8', 'AFF8h', 'F8'});



% Specify the following parameters
NUM_DIP_WANTED = 10000;
NUM_SOURCES = 1;

num_dip = size(lf.leadfield, 2);
rng(0, 'twister');
%random_columns = randperm(num_dip, NUM_DIP_WANTED);

random_columns =  lf_get_source_spaced(lf, NUM_DIP_WANTED, 3);
L_r = lf.leadfield(:, random_columns);

for i=1:NUM_DIP_WANTED
    L_r(:,i) = squeeze(lf.leadfield(:,random_columns(i),:)) * lf.orientation(random_columns(i),:)';
end

remaining_columns = setdiff(1:num_dip, random_columns);
L_o = lf.leadfield(:, remaining_columns);
W_s = weighting(L_r, 64);
TOT_SOURCES = 1;

W1 = zeros(TOT_SOURCES, 1);
W2 = zeros(TOT_SOURCES, 1);
EMD = zeros(TOT_SOURCES, 1);

true_sources_pos = zeros(TOT_SOURCES, 1);
inv_sources_pos = zeros(TOT_SOURCES, 1);
inv_sources = cell(TOT_SOURCES, 1);

active_dipole_1 = cell(TOT_SOURCES, 1); 
active_dipole_2 = cell(TOT_SOURCES, 1); 
value = 1e-6; 


for k = 1:TOT_SOURCES
    disp(['Simulation:', num2str(k)])
    source_vec = zeros(size(L_o,2), 1);
    source_location = find_source_location(L_o);
    source_vec(source_location) = 1;
    true_sources_pos(k) = source_location;

    data = L_o * source_vec;
    sigma = 0.05;
    alpha = .05;
    noise_E = sigma.^2 .* eye(length(data)) * randn(size(data)); % Covariance matrix of measurement noise
    noisy_data = data + noise_E;
  
    W_s_inv = diag(1./diag(W_s));
    fidelity_error = 2;
    while fidelity_error > 1.3
        [f_inv, status] = l1_ls(L_r * W_s_inv, noisy_data, alpha, 1e-3);
        fs_inv = W_s_inv * f_inv;
        fidelity_error = (norm(L_r*(fs_inv) - noisy_data)) / norm(noise_E);
        alpha = 0.5*alpha
    end
    %[f_inv, status] = l1_ls(L_r * W_s_inv, data, 1e-3, 1e-4);
    %fs_inv = W_s_inv * f_inv;
    [max_fs_inv, max_index] = max(fs_inv);
    inv_sources_pos(k) = max_index;
    ts_50(k) = length(find(fs_inv > 0.5 * max(fs_inv)));
    ts_25(k) = length(find(fs_inv > 0.25 * max(fs_inv)));
  

    
    F1 = lf.pos(remaining_columns(true_sources_pos(k)), :); 
    F2 = lf.pos(random_columns(inv_sources_pos(k)), :);  
    W1 = 1;  
    W2 = fs_inv(max_index);  

    W1 = W1 / sum(W1);
    W2 = W2 / sum(W2);


    [~, fval] = emd(F1, F2, W1, W2, @gdf);
    EMD(k) = fval;

 
end




fprintf('||noise||: %f\n', norm(noise_E));
fprintf('||L_r*(fs_inv) - noisy_data)||: %f\n', norm(L_r*(fs_inv) - noisy_data));

EEG_coordinates = [lf.chanlocs.X; lf.chanlocs.Y; lf.chanlocs.Z]';

percentage = (norm(data - noisy_data)/ norm(data))*100;

fprintf('Noise Percentage: %f\n', percentage);


fprintf('Earth Mover''s Distance: %f\n', mean(EMD));


euclidean_distances = sqrt(sum((F2 - F1).^2, 2));
min_direct_distance = min(euclidean_distances);


fprintf('Minimum Direct Euclidean Distance: %f\n', min_direct_distance);
fprintf('EMD: %f\n', mean(EMD));

fprintf('Mean ts_50: %f\n', mean(ts_50));
fprintf('Mean ts_25: %f\n', mean(ts_25));
