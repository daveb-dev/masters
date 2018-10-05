close all
filesNums = {'01','02','05','06','09','12'};
rat = struct('tslices',[],'area',[],'skulls',{},'initial',[],'time',[],...
    'meancells',[],'meanarea',[]);
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    days = size(cells,4);
    rat(n).tslices = zeros(days,16);
    rat(n).area = zeros(days,16);
    rat(n).skulls = cell(1,16);
    
    for z = 1:16
        % Slices with tumor cells
        for t = 1:days
            c = sum(sum(cells(:,:,z,t)));
            a = sum(sum(cells(:,:,z,t)>0));
            rat(n).tslices(t,z) = c;
            rat(n).area(t,z) = a;
            
            % Time
            rat(n).time(t) = time(t);
        end
        
        % Skull coordinates for eachs slice
        [row,col,~] = find(skull(:,:,z) == 1);
        rat(n).skulls{z} = [row col];

    end
    
    % Number of tumor cells vs time
    figure(1)
    sumSlices = sum(rat(n).tslices,2);
    plot(sumSlices); hold on
        
    % "Volume" (# of pixels) of tumor cells vs time
    figure(2)
    sumArea = sum(rat(n).area,2);
    plot(sumArea); hold on
    
    % Calculating mean tumor cells
    if size(sumSlices,1) > size(rat(1).meancells,1)
        rat(1).meancells(size(sumSlices,1),n) = 0;
        rat(1).meanarea(size(sumSlices,1),n) = 0;
    elseif size(sumSlices,1) < size(rat(1).meancells,1)
        sumSlices(size(rat(1).meancells,1),1) = 0;
        sumArea(size(rat(1).meancells,1),1) = 0;
    end
    rat(1).meancells = [rat(1).meancells sumSlices];
    rat(1).meanarea = [rat(1).meanarea sumArea];
    
    % Plot evolution of tumor in brain
    f = gcf;
    [~,maxlayer] = max(sum(rat(n).tslices));  % Slice with most cells
    rat(n).initial = cells(:,:,maxlayer,1);   % Initial value
    if days <= 6
        b = 3;
    else
        b = 4;
    end
    tumor = cells(:,:,maxlayer,:);
    mincell = min(tumor(tumor>0));  maxcell = max(max(max(max(tumor))));
    for t = 1:days
        
        % Turn data into RGB for layering on brain
        tumor = cells(:,:,maxlayer,t);
        [B,~] = real2rgb(tumor,'jet',[mincell maxcell]);
        B3 = B(:,:,3); B3(tumor==0) = 0; B(:,:,3) = B3;
        
        % Turn off brain data where tumor goes, convert to RGB
        brain = anatomical(:,:,maxlayer,t);
        [C,~] = real2rgb(brain,'bone',[0 45]);
        t3 = [tumor tumor tumor]; t3 = reshape(t3,41,61,3);
        C(t3 > 0 ) = 0;
        
        % Add images together and plot
        D = B+C;
        figure(f.Number+n)
        subplot(2,b,t)
        imshow(D); hold on
        title(strcat('Day',num2str(time(t))));
        ax = gca;
        ax.Position = ax.Position + [-.03 0 0 0]*(b-mod(t-1,b));
        
        % Skull outline
        scatter(rat(n).skulls{maxlayer}(:,2),rat(n).skulls{maxlayer}(:,1),'k.')
        
    end
    set(gcf,'NextPlot','add'); axes;
    h = title('Example of Tumor Growth');
    cbar = colorbar('eastoutside');
    caxis([mincell maxcell]); colormap jet;
    set(gca,'Visible','off');
    set(h,'Visible','on');
    pos=get(cbar,'Position');
    set(cbar,'Position',pos+[0.06,0,0,0]);
    
end

figure(1)
rat(1).meancells(rat(1).meancells == 0) = nan;
plot(mean(rat(1).meancells,2,'omitnan'),'Linewidth',2);
title('Evoluation of Total Tumor Cells in Rat Brain')
xlabel('time (days)')
ylabel('Number of Tumor Cells')
legend([filesNums 'mean'],'Location','Northwest');

figure(2)
rat(1).meanarea(rat(1).meanarea == 0) = nan;
plot(mean(rat(1).meanarea,2,'omitnan'),'Linewidth',2);
title('Evoluation of Total Tumor Volume in Rat Brain')
xlabel('time (days)')
ylabel('Volume of Tumor Cells (pixels with cells)')
legend([filesNums 'mean'],'Location','Northwest');

save('finaldata.mat','rat');
