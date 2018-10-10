close all
filesNums = {'01','02','05','06','09','12'};
rat = struct('tslices',[],'area',[],'initial',[],'time',[],...
    'sbound',{},'tbound',{},'meanarea',[],'maxlayer',[]);

%% Number of cells in slice z at time t
close all
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    days = size(cells,4);
    rat(n).tslices = zeros(days,16);
    for z = 1:16
        for t = 1:days
            c = sum(sum(cells(:,:,z,t)));
            rat(n).tslices(t,z) = c;
        end
    end
    % Number of tumor cells vs time
    sumSlices = sum(rat(n).tslices,2);
    plot(time,sumSlices,'o-','Linewidth',2); hold on
end
xlabel('time (days)','FontSize',14)
ylabel('Number of Tumor Cells','FontSize',14)
legend(filesNums,'Location','Northwest','FontSize',14);
set(get(legend,'Title'),'String','Rat ID')
saveas(gcf,'images/numtumorcells','png');
%% Time
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    rat(n).time = time;
end
%% Initial condition based on max initial count
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    [~,rat(n).maxlayer] = max(rat(n).tslices(1,:));  % Slice with most cells
    rat(n).initial = cells(:,:,rat(n).maxlayer,1);   % Initial value
end
%% Area
close all
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    days = size(cells,4);
    rat(n).area = zeros(days,16);
    for z = 1:16
        for t = 1:days
            a = sum(sum(cells(:,:,z,t)>0));
            rat(n).area(t,z) = a;
        end
    end
    area = rat(n).area(:,rat(n).maxlayer);
    plot(rat(n).time,area,'o-','Linewidth',2); hold on
end
xlabel('time (days)','FontSize',14)
ylabel('# of Voxels with Tumor Cells)','FontSize',14)
legend(filesNums,'Location','Northwest','FontSize',14);
set(get(legend,'Title'),'String','Rat ID')
saveas(gcf,'images/areatumorcells','png');
%% Volume
close all
for n = 1:length(filesNums)
    volume = sum(rat(n).area,2);
    plot(rat(n).time,volume,'o-','Linewidth',2); hold on
end
xlabel('time (days)','FontSize',14)
ylabel('# of Voxels with Tumor Cells)','FontSize',14)
legend(filesNums,'Location','Northwest','FontSize',14);
set(get(legend,'Title'),'String','Rat ID')
saveas(gcf,'images/voxtumorcells','png');
%% Boundaries
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    [~,z] = max(rat(n).tslices(1,:));
    %rat(n).skulls = cell(1,16);
    rat(n).sbound = cell(1,length(rat(n).time));
    rat(n).tbound = cell(1,length(rat(n).time));
    for t = 1:length(rat(n).time)
        [row,col,~] = find(skull(:,:,z) == 1);
        k = boundary(row,col);
        rat(n).sbound{t} = [col(k) 41-row(k)];
        [row,col,~] = find(cells(:,:,z,t) > 0);
        k = boundary(row,col,1);
        rat(n).tbound{t} = [col(k) 41-row(k)];
    end
end
%% Montages
for n = 1:length(filesNums)
    close all
    load(strcat('W',filesNums{n},'_model_data.mat'));
    days = size(cells,4);
    z = rat(n).maxlayer;
    tumor = cells(:,:,z,:); % Slice with most cells
    mincell = min(tumor(tumor>0));  maxcell = max(max(max(max(tumor))));
    imds = [];
    for t = 1:days
        % Turn data into RGB for layering on brain
        tumor = cells(:,:,z,t);
        [B,~] = real2rgb(tumor,'jet',[mincell maxcell]);
        B3 = B(:,:,3); B3(tumor==0) = 0; B(:,:,3) = B3;
        
        % Turn off brain data where tumor goes, convert to RGB
        brain = anatomical(:,:,z,t);
        [C,~] = real2rgb(brain,'bone',[0 45]);
        t3 = [tumor tumor tumor]; t3 = reshape(t3,41,61,3);
        C(t3 > 0 ) = 0;
        
        % Skull and tumor outline
        xy = []; 
        s = rat(n).sbound{t};
        for i = 1:length(s)-1
            xy = [xy(1:end-1,:); bresenham(s(i,1),s(i,2),s(i+1,1),s(i+1,2))];
        end
        xy(end+1,:) = xy(end,:);
        tu = rat(n).tbound{t};
        for i = 1:length(tu)-1
            xy = [xy(1:end-1,:); bresenham(tu(i,1),tu(i,2),tu(i+1,1),tu(i+1,2))];
        end
        sout = (41-xy(:,2))+41*(xy(:,1)-1);
        sout = [sout; sout+41*61; sout+41*61*2];
        
        % Add images together and cat for montage
        D = B+C;
        D(sout) = 0;
        imds = cat(4,imds,D);
    end
    montage(imds,'Size',[NaN 3]); hold on
    cbar = colorbar('eastoutside'); hold on
    ylabel(cbar,'Number of Tumor Cells in Pixel','FontSize',14)
    caxis([mincell maxcell]); colormap jet;
    saveas(gcf,strcat('images/Montage',filesNums{n}),'png');
end
%% Test Mean

%% Last
% for n = 1:length(filesNums)
%     load(strcat('W',filesNums{n},'_model_data.mat'));
%     days = size(cells,4);
%     rat(n).tslices = zeros(days,16);
%     rat(n).area = zeros(days,16);
%     rat(n).skulls = cell(1,16);
%     rat(n).time = time;
%     for z = 1:16
%         % Slices with tumor cells
%         for t = 1:days
%             c = sum(sum(cells(:,:,z,t)));
%             rat(n).tslices(t,z) = c;
%             a = sum(sum(cells(:,:,z,t)>0));
%             rat(n).area(t,z) = a;
%         end
%         
%         % Skull & tumor boundary for each slice
%         [row,col,~] = find(skull(:,:,z) == 1);
%         k = boundary(row,col);
%         rat(n).skulls{z} = [row(k) col(k)];
%         [row,col,~] = find(cells(:,:,z) == 1);
%         k = boundary(row,col);
%         rat(n).skulls{z} = [row(k) col(k)];
%     end
%     rat(n).sbound = rat(n).skulls{maxlayer};
%     
%     % Number of tumor cells vs time
%     figure(1)
%     sumSlices = sum(rat(n).tslices,2);
%     plot(time,sumSlices,'o-'); hold on
%         
%     % "Volume" (# of pixels) of tumor cells vs time
%     figure(2)
%     sumArea = sum(rat(n).area,2);
%     plot(time,sumArea,'o-'); hold on
%     
%     % Plot evolution of tumor in brain
%     [~,maxlayer] = max(rat(n).tslices(1,:));  % Slice with most cells
%     rat(n).initial = cells(:,:,maxlayer,1);   % Initial value
%     if days <= 6
%         b = 3;
%     else
%         b = 4;
%     end
%     tumor = cells(:,:,maxlayer,:);
%     mincell = min(tumor(tumor>0));  maxcell = max(max(max(max(tumor))));
%     imds = [];
%     for t = 1:days
%         % Turn data into RGB for layering on brain
%         tumor = cells(:,:,maxlayer,t);
%         [B,~] = real2rgb(tumor,'jet',[mincell maxcell]);
%         B3 = B(:,:,3); B3(tumor==0) = 0; B(:,:,3) = B3;
%         
%         % Turn off brain data where tumor goes, convert to RGB
%         brain = anatomical(:,:,maxlayer,t);
%         [C,~] = real2rgb(brain,'bone',[0 45]);
%         t3 = [tumor tumor tumor]; t3 = reshape(t3,41,61,3);
%         C(t3 > 0 ) = 0;
%         
%         % Add images together and cat for montage
%         D = B+C;
%         imds = cat(4,imds,D);
%     end
%     f = gcf;
%     figure(f.Number+n)
%     montage(imds); hold on
%     cbar = colorbar('eastoutside'); hold on
%     caxis([mincell maxcell]); colormap jet;
%     
% end
% 
% figure(1)
% rat(1).meancells(rat(1).meancells == 0) = nan;
% plot(mean(rat(1).meancells,2,'omitnan'),'Linewidth',2);
% title('Evoluation of Total Tumor Cells in Rat Brain')
% xlabel('time (days)')
% ylabel('Number of Tumor Cells')
% legend([filesNums 'mean'],'Location','Northwest');
% 
% figure(2)
% rat(1).meanarea(rat(1).meanarea == 0) = nan;
% plot(mean(rat(1).meanarea,2,'omitnan'),'Linewidth',2);
% title('Evoluation of Total Tumor Volume in Rat Brain')
% xlabel('time (days)')
% ylabel('Volume of Tumor Cells (pixels with cells)')
% legend([filesNums 'mean'],'Location','Northwest');

save('finaldata.mat','rat');
