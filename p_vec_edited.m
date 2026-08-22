clear all; clc; 
lf = lf_generate_fromnyhead('montage', 'S64');
%lf = lf_generate_fromnyhead();
% lf = lf_generate_fromnyhead('labels', {'Nz', 'LPA', 'RPA', 'PPO9h', 'PO7', 'PO9', 'POO9h', ...
 %         'O1', 'I1', 'OI1h', 'Oz', 'Iz', 'OI2h', 'O2', 'I2', 'POO10h', ...
  %          'PO8', 'PO10', 'PPO10h', 'PPO7h', 'PO5h', 'PPO3h', 'POO5', ...
   %         'PO3h', 'PPO1h', 'POO3', 'Pz', 'POz', 'PPO2h', 'POO4',...
    %        'PO4h', 'PPO4h', 'POO6', 'PO6h', 'PPO8h', 'FFT9h', 'FC5', ...
     %       'FT7', 'FT9', 'FCC5h', 'FTT7h', 'C5', 'T7', 'CCP5h', ...
      %      'TTP7h', 'CP5', 'TP7', 'TPP7h', 'P5', 'P7', 'P9', 'P6',...
       %     'P8', 'P10', 'TPP8h', 'CP6', 'TP8', 'CCP6h', 'TTP8h',...
        %    'C6', 'T8', 'FCC6h', 'FTT8h', 'FC6', 'FT8', 'FT10', ...
         %   'FFT10h', 'FC1', 'FC3', 'FCC1h', 'FCC3h', 'C1', 'C3', ...
          %  'CCP1h', 'CCP3h', 'CPz', 'CP1', 'CP3', 'CPP1h', 'CPP3h',...
           % 'CPP5h', 'P1', 'P3', 'P2', 'P4', 'CPP2h', 'CPP4h', 'CPP6h',...
%            'CP2', 'CP4', 'CCP2h', 'CCP4h', 'Cz', 'C2', 'C4', 'FCC2h', ...
 %           'FCC4h', 'FC2', 'FC4', 'FFT7h', 'FFC5h', 'F5', 'FFC3h', ...
  %          'F3', 'FFC1h', 'F1', 'FCz', 'Fz', 'FFC2h', 'F2', 'FFC4h',...
   %         'F4', 'FFC6h', 'F6', 'FFT8h', 'F7', 'AFF7h', 'AF5h', 'AF7', ...
    %        'AF3h', 'Fp1', 'AFF1h', 'AFF2h', 'AFz', 'Fpz', 'AF4h', 'Fp2',...
     %       'AF6h', 'AF8', 'AFF8h', 'F8'});

% ORIENTED DIPOLE ONE SOURCE WITHOUT NOISE
rng(0, 'twister');
random_numbers = rand(4, 1);  
random_numbers = random_numbers / norm(random_numbers);  
d1 = random_numbers(1);
d2 = random_numbers(2);
d3 = random_numbers(3);
d4 = random_numbers(4);


N = size(lf.orientation, 1);
alpha = pi/6;
c1 = cos(alpha);
c2 = sin(alpha);
v1 = zeros(N, 3);
v2 = zeros(N, 3);
v3 = zeros(N, 3);
v4 = zeros(N, 3);

for i = 1:N
    orientation_vector = lf.orientation(i, :) / norm(lf.orientation(i, :));
    a_vector = [-orientation_vector(2), orientation_vector(1), 0] / norm([-orientation_vector(2), orientation_vector(1), 0]);
    perpendicular_vector = cross(orientation_vector, a_vector) / norm(cross(orientation_vector, a_vector)); 

    v1(i,:) = c1 * orientation_vector + c2 * a_vector;
    v2(i,:) = c1 * orientation_vector - c2 * a_vector;
    v3(i,:) = c1 * orientation_vector + c2 * perpendicular_vector;
    v4(i,:) = c1 * orientation_vector - c2 * perpendicular_vector;
end


num_eeg = size(lf.leadfield, 1);
num_dip = size(lf.leadfield, 2);
L_v1 = zeros(num_eeg, num_dip);
L_v2 = zeros(num_eeg, num_dip);
L_v3 = zeros(num_eeg, num_dip);
L_v4 = zeros(num_eeg, num_dip);

for j = 1:num_dip
    L = squeeze(lf.leadfield(:, j, :));
    L_v1(:, j) = L * v1(j, :)';
    L_v2(:, j) = L * v2(j, :)';
    L_v3(:, j) = L * v3(j, :)';
    L_v4(:, j) = L * v4(j, :)';
end

L_new = [L_v1, L_v2, L_v3, L_v4];


NUM_DIP_WANTED = 2500; 
random_columns = lf_get_source_spaced(lf, NUM_DIP_WANTED, 3);
Linv = [L_v1(:, random_columns), L_v2(:, random_columns), L_v3(:, random_columns), L_v4(:, random_columns)];
W_s = weighting(Linv, 64);
TOT_SOURCES = 1;

W1 = zeros(TOT_SOURCES, 1);
W2 = zeros(TOT_SOURCES, 1);
EMD = zeros(TOT_SOURCES, 1);
active_dipole_1 = cell(TOT_SOURCES, 1); 
active_dipole_2 = cell(TOT_SOURCES, 1);

true_sources_pos = zeros(TOT_SOURCES, 1);
inv_sources_pos = zeros(TOT_SOURCES, 1);
all_inv_sources_pos = cell(TOT_SOURCES, 1);
W_s_inv = diag(1./diag(W_s));

for k = 1:TOT_SOURCES
    disp(['Simulation:', num2str(k)])
    source_vec = zeros(size(L_new, 2), 1);
    source_location = find_source_location(L_v1);    
    source_vec(source_location) = d1;
    source_vec(source_location + num_dip) = d2;
    source_vec(source_location + 2*num_dip) = d3;
    source_vec(source_location + 3*num_dip) = d4;
    true_sources_pos(k) = source_location;

    
    data = L_new * source_vec;
    sigma = 0.04;
    noise_E = sigma.^2 .* eye(length(data)) * randn(size(data)); 
    noisy_data = data + noise_E;
    beta= 0.025;  
    fidelity_error = 2;  

    while fidelity_error > 1.3 
    prev_error = fidelity_error;  
    [f_inv, status] = l1_ls_nonneg(Linv * W_s_inv, noisy_data, beta, 1e-3);
    fs_inv = W_s_inv * f_inv;
        
    f_inverse = reshape(fs_inv,[2500,4]);
    f_norm = sum(sqrt(f_inverse.^2),2);
    fidelity_error = (norm(Linv*(fs_inv) - noisy_data)) / norm(noise_E);
    beta = 0.9*beta
    end

 

    
    %[f_inv, status] = l1_ls_nonneg(Linv * W_s_inv, data, 1e-3, 1e-4);
    %fs_inv = W_s_inv * f_inv;
    %f_inverse = reshape(fs_inv, [NUM_DIP_WANTED, 4]);
    %f_norm = sum(sqrt(f_inverse.^2), 2);
    [max_fs_inv, max_index] = max(f_norm); 
   
    inv_sources_pos(k) = max_index;

    F1 = lf.pos(true_sources_pos(k), :);  
    F2 = lf.pos(random_columns, : );  
    W1 = 1;  
    W2 = f_norm; 

    W1 = W1 / sum(W1);
    W2 = W2 / sum(W2);

  
    [~, fval] = emd(F1, F2, W1, W2, @gdf);
    EMD(k) = fval;
    
    ts_50(k) = length(find(f_norm > .5*f_norm(max_index)));
    ts_25(k) = length(find(f_norm > .25*f_norm(max_index)));    
    threshold_1 = 0.25;
    threshold_2 = 0.5;
    close_value_1 = max_fs_inv * threshold_1;
    close_value_2 = max_fs_inv * threshold_2;
    close_indices_1 = find(fs_inv>= close_value_1);
    close_indices_2 = find(fs_inv >= close_value_2);
    active_dipole_1{k} = close_indices_1;
    active_dipole_2{k} = close_indices_2;


end

total_active_dipoles_1 = 0;
total_active_dipoles_2 = 0;

for k = 1:TOT_SOURCES
    num_active_dipoles_1 = length(active_dipole_1{k});
    total_active_dipoles_1 = total_active_dipoles_1 + num_active_dipoles_1;
end

for k = 1:TOT_SOURCES
    num_active_dipoles_2 = length(active_dipole_2{k});
    total_active_dipoles_2 = total_active_dipoles_2 + num_active_dipoles_2;
end



spatial_disp_25 = total_active_dipoles_1 / TOT_SOURCES;
fprintf('Spacial Dispersion with 0.25: %f\n', spatial_disp_25);

spatial_disp_50 = total_active_dipoles_2 / TOT_SOURCES;
fprintf('Spacial Dispersion with 0.50: %f\n', spatial_disp_50);


% Print results
fprintf('Mean Earth Mover''s Distance: %f\n', mean(EMD));
fprintf('Noise Percentage: %f\n', (norm(data - noisy_data)/ norm(data)) * 100);
fprintf('||noise||: %f\n', norm(noise_E));
fprintf('||Linv*(fs_inv) - noisy_data)||: %f\n', norm(Linv*(fs_inv) - noisy_data));
fprintf('Mean ts_50: %f\n', mean(ts_50));
fprintf('Mean ts_25: %f\n', mean(ts_25));