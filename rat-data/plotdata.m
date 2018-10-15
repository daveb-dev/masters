close all
filesNums = {'01','02','05','06','09','12'};
rat = struct('tslices',[],'area',[],'initial',[],'time',[],...
    'sbound',[],'tbound',{},'meanarea',[],'maxlayer',[]);

% Number of cells in slice z at time t
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

% Time
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    rat(n).time = time;
end

% Initial condition & data based on max initial count
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    [~,rat(n).maxlayer] = max(rat(n).tslices(1,:));  % Slice with most cells
    ic = cells(:,:,rat(n).maxlayer,1);   % Initial value rat(n).initial
    fstr = strcat('rat',filesNums{n});
    eval(strcat('save(''./',fstr,'/ic.mat'',''ic'');'));
    for t = 2:length(rat(n).time)
        tumor = cells(:,:,rat(n).maxlayer,t);  
        days = num2str(rat(n).time(t)-rat(n).time(1));
        eval(strcat('save(''./',fstr,'/tumor_t',days,'.mat'',''ic'');'));
    end
end

% Boundaries
for n = 1:length(filesNums)
    load(strcat('W',filesNums{n},'_model_data.mat'));
    [~,z] = max(rat(n).tslices(1,:));
    
    % Skull
    skull_out = [];
    [row,col,~] = find(skull(:,:,z) == 1);
    k = boundary(row,col);
    s = [col(k) 41-row(k)];
    for i = 1:length(s)-1
        skull_out = [skull_out(1:end-1,:); bresenham(s(i,1),s(i,2),s(i+1,1),s(i+1,2))];
    end
    rat(n).sbound = skull_out;
    fstr = strcat('rat',filesNums{n});
    eval(strcat('save(''./',fstr,'/skull_out.mat'',''skull_out'');'));
    
    % Tumor
    tumor_out = [];
    [row,col,~] = find(cells(:,:,z,1) > 0);
    k = boundary(row,col,1);
    tu = [col(k) 41-row(k)];
    for i = 1:length(tu)-1
        tumor_out = [tumor_out(1:end-1,:); bresenham(tu(i,1),tu(i,2),tu(i+1,1),tu(i+1,2))];
    end
    rat(n).tbound = tumor_out;
    eval(strcat('save(''./',fstr,'/tumor_out.mat'',''tumor_out'');'));
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
        xy = [rat(n).sbound{t}; rat(n).tbound{t}];
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


save('finaldata.mat','rat');
