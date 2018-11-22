 
% calc_echo_spacing_philips_mod.m
%
% Formula from http://www.spinozacentre.nl/wiki/index.php/NeuroWiki:Current_developments#B0_correction
% modified after suggestions from Pieter Buur and Nikolaus Weiskopf; also thanks to Lawrie MacKay.
% Please check results carefully, if possible with physicist of the Philips scanner where data were acquired from.
% 
% For Philips Achieva, the best equation for calculating "effective" echo spacing 
% (corrected for SENSE factor) in milliseconds would seem to be:
%  
% echo spacing in msec = (1000 * wfs)/(434.215 * (etl)) (etl = EPI factor
% + 1)
% 
% We have:
% - matrix size in the phase encoding direction 
% - EPI factor (EPI factor + 1 = echo train length (ETL))
% - water-fat shift per pixel
% - field strength in Tesla

% % We calculate:
% - total bandwidth
% - water-fat shift in Hz
% - bandwidth per pixel in Hz
% - echo spacing in sec
% - echo spacing in msec

% We assume:
% - that the SENSE factor is incorporated in Philips' EPI factor
% - that the water-fat difference in parts-per-million (ppm) is: 3.35 (from Haacke)
% 
% Use for own risk. Please don't use for clinical applications. 
% Hester Breman, 23-08-2013; latest update: 19-11-14
 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variables, please adapt for each dataset (can be found in *.PAR file) ~~~~~~~~

matrixsize_phase_enc_dir    = 80; % matrix size in y direction (if that is phase encoding direction)
epifactor 					= 43; % EPI factor
water_fat_shift_pixel 		= 12.0; % water-fat shift per pixel
fieldstrength_tesla 		= 3.0;  % magnetic field strength (T)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ not change below line ~~~~~~~~

water_fat_diff_ppm          = 3.35;  
resonance_freq_mhz_tesla    = 42.576; % gyromagnetic ratio for proton (1H)
echo_train_length           = epifactor + 1 
water_fat_shift_hz          = fieldstrength_tesla * water_fat_diff_ppm * resonance_freq_mhz_tesla % water_fat_shift_hz 3T = 434.215 Hz
BW_hz_pixel                 = water_fat_shift_hz / water_fat_shift_pixel
totBW                       = BW_hz_pixel * echo_train_length
echo_spacing_sec            = 1/totBW 
echo_spacing_msec           = echo_spacing_sec * 1000

