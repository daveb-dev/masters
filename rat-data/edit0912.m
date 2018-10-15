% This code was used to edit rats 09 and 12 data so that imaging values
% were comparable to other images

clear
load('W09_model_data.mat');
anatomical = anatomical_original;
divisor = max(max(max(anatomical(:,:,12,:))))/65;
anatomical = anatomical/divisor;
anatomical(:,:,:,6) = anatomical(:,:,:,6)*2;
anatomical(:,:,:,2:end) = anatomical(:,:,:,2:end)*1.2;
save('W09_model_data.mat');


clear
load('W12_model_data.mat');
anatomical = anatomical_original;
anatomical(:,:,:,3) = anatomical(:,:,:,4);
divisor = max(max(max(anatomical(:,:,10,:))))/100;
anatomical = anatomical/divisor;
save('W12_model_data.mat');

clear