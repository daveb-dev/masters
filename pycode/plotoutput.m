% close all
% filesNums = {'01','02','05','06','09','12'};
% rat = struct('tslices',[],'area',[],'initial',[],'time',[],...
%     'sbound',[],'tbound',{},'meanarea',[],'maxlayer',[]);

%% initial
close all
imds = {imread('./output/rat05_day0.png')};
montage(imds,'Size',[NaN 1]); hold on
cbar = colorbar('eastoutside','FontSize',15); hold on
ylabel(cbar,'Fraction of Cancerous Cells','FontSize',20)
caxis([0 1]); colormap jet;
saveas(gcf,'./output/rat05_day0_2','png');

%% Montages of tumor growth
for n = 1
    close all
    days = [0 2 4 5 6 9]; %size(cells,4);
    imds = {};
    for t = 1:length(days)
        imds{t} = imread(strcat('./output/lehe_day',num2str(days(t)),'.png'));
    end
    montage(imds,'Size',[NaN 3]); hold on
    cbar = colorbar('eastoutside','FontSize',15); hold on
    ylabel(cbar,'Fraction of Cancerous Cells','FontSize',20)
    caxis([0 1]); colormap jet;
    saveas(gcf,'./output/Montage_lehe','png');
    
%     imds = {};
%     for t = 1:length(days)
%         imds{t} = imread(strcat('./output/le_day',num2str(days(t)),'.png'));
%     end
%     montage(imds,'Size',[NaN 3]); hold on
%     cbar = colorbar('eastoutside','FontSize',15); hold on
%     ylabel(cbar,'Fraction of Cancerous Cells','FontSize',20)
%     caxis([0 1]); colormap jet;
%     saveas(gcf,'./output/Montage_le','png');
%     
%     imds = {};
%     for t = 1:length(days)
%         imds{t} = imread(strcat('./output/he_day',num2str(days(t)),'.png'));
%     end
%     montage(imds,'Size',[NaN 3]); hold on
%     cbar = colorbar('eastoutside','FontSize',15); hold on
%     ylabel(cbar,'Fraction of Cancerous Cells','FontSize',20)
%     caxis([0 1]); colormap jet;
%     saveas(gcf,'./output/Montage_he','png');
end


%% K/D fields
close all
imds = {imread('./output/le_D0.png')};
montage(imds,'Size',[NaN 1]); hold on
cbar = colorbar('eastoutside','FontSize',15); hold on
ylabel(cbar,'D_0, diffusion coefficient field','FontSize',20)
caxis([0 5]); colormap jet;
saveas(gcf,'./output/le_D0_map','png');

close all
imds = {imread('./output/he_D0.png')};
montage(imds,'Size',[NaN 1]); hold on
cbar = colorbar('eastoutside','FontSize',15); hold on
ylabel(cbar,'Fraction of Cancerous Cells','FontSize',20)
caxis([0 5]); colormap jet;
saveas(gcf,'./output/he_D0_map','png');

close all
imds = {imread('./output/le_k0.png')};
montage(imds,'Size',[NaN 1]); hold on
cbar = colorbar('eastoutside','FontSize',15); hold on
ylabel(cbar,'k_0, growth coefficient field','FontSize',20)
caxis([0 10]); colormap jet;
saveas(gcf,'./output/le_k0_map','png');

close all
imds = {imread('./output/he_k0.png')};
montage(imds,'Size',[NaN 1]); hold on
cbar = colorbar('eastoutside','FontSize',15); hold on
ylabel(cbar,'k_0, growth coefficient field','FontSize',20)
caxis([0 10]); colormap jet;
saveas(gcf,'./output/he_k0_map','png');

%%
close all
days = [0 2 4 5 6 9]; %size(cells,4);
true_num_cells = zeros(1,length(days));
le_num_cells = zeros(1,length(days));
he_num_cells = zeros(1,length(days));
true_num_cells2 = zeros(1,length(days));

true_area_cells = zeros(1,length(days));
le_area_cells = zeros(1,length(days));
he_area_cells = zeros(1,length(days));

lim = .01;

for t = 1:6
    data = h5read('./output/rat05le/nosteps.h5',...
        strcat('/VisualisationVector/',num2str(2*t-1)));
    true_num_cells(t) = sum(data)*50970;
    true_area_cells(t) = sum(data>.01);
    
    data = h5read('./output/rat05le/nosteps.h5',...
        strcat('/VisualisationVector/',num2str(2*(t-1))));
    le_num_cells(t) = sum(data)*50970;
    le_area_cells(t) = sum(data>lim);
    
    data = h5read('./output/rat05he/nosteps.h5',...
        strcat('/VisualisationVector/',num2str(2*t-1)));
    true_num_cells2(t) = sum(data)*50970;
    
    data = h5read('./output/rat05he/nosteps.h5',...
        strcat('/VisualisationVector/',num2str(2*(t-1))));
    he_num_cells(t) = sum(data)*50970;
    he_area_cells(t) = sum(data>lim);
    
end

figure(1)
plot(days,true_num_cells,'o',days,le_num_cells,'x',days,he_num_cells,...
    '^','LineWidth',2);
legend('Actual Data','LE Model','HE Model','Location','Northwest');
ylabel('Number of tumor cells');
xlabel('Day');
xlim([0 10])

saveas(gcf,'./output/num_cells','png');

figure(2)
plot(days,true_area_cells,'o',days,le_area_cells,'x',days,he_area_cells,...
    '^','LineWidth',2);
legend('Actual Data','LE Model','HE Model','Location','Northwest');
ylabel('"Area" of tumor cells');
xlabel('Day');
xlim([0 10])

saveas(gcf,'./output/area_cells','png');
close all

