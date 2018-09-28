clear all; close all
cd /h2/xoab/repos/research/code/i_hooke_variance_ml/outputData
filtChain_ml
rawChain_ml
filt = ip_ml_last_filtChain_unified;
filt = filt(:,1);
C = who('*rawChain_unified');
for c = 1:length(C)
    D = strsplit(C{c},'_');
     Ch(c) = str2num(D{3});
end
lastch = strcat('ip_ml_',num2str(max(Ch)),'_rawChain_unified');
raw = eval(lastch);
raw = raw(:,1);
cd ..

blue = [0 .447 .741];
% % Density plots
band = 10000;
band2 = 20000;
pts = linspace(1000,300000,1000);
[f,xi] = ksdensity(filt,pts,'function','pdf','Bandwidth',band);
plot(xi,f,'-r','linewidth',3); hold on
x = [1000 300000];
Epri = pdf('unif',x,x(1),x(2));
plot(x,Epri,'Color',blue,'Linewidth',3);
ylim([0 1.25*max(f)]);
xlabel('E (MPa)');
ylabel('KDE');
title('Posterior for E');
legend('Posterior Distribution','Prior Distribution','Location','NorthEast');
print -dpng le_E_ml.png
pause; clf;

% [f,xi] = ksdensity(raw,pts,'function','pdf','Bandwidth',band);
% plot(xi,f,'-r','linewidth',3); hold on
% plot(x,Epri,'Color',blue,'Linewidth',3);
% ylim([0 max(f)*1.25])
% xlabel('E (MPa)');
% ylabel('KDE');
% title('Posterior for E');
% legend('Posterior Distribution','Prior Distribution','Location','NorthEast');
% pause; clf;

% % Chain plots
% plot(raw); hold on
% plot(filt,'r')
% ylabel('\theta=E','fontname', 'Times', 'fontsize',12);
% xlabel('Number of positions','fontname', 'Times', 'fontsize',12);
% title('DRAM MCMC Chain Positions (filtered)','fontname', 'Times', 'fontsize',12);
% legend('Raw','Filtered')
% print -dpng E_chain_pos_filt_ml.png
% pause; clf; 
% 
% 
% % % Autocorrelation plots 
% % % RAW and FILTERED
% [ACFr,lagsr,~] = autocorr(raw); hold on
% [ACFf,lagsf,~] = autocorr(filt);
% plot(lagsr,ACFr,'o-',lagsf,ACFf,'ro-','Linewidth',2);
% ylabel('Autocorrelation for \theta=E','fontname', 'Times', 'fontsize',12);
% xlabel('Lag');
% title('Parameter Autocorrelation');
% h=legend('Raw chain','Filtered chain','location','northeast');
% print -dpng E_autocorr_ml.png
% pause; close 
% 
% % 
% % % Covariance matrix 
% cov_matrix = cov(ip_ml_last_filtChain_unified)

% CDF
[f,xi] = ksdensity(filt,pts,'function','cdf','Bandwidth',band2);
plot(xi,f,'-r','linewidth',3); hold on
xlabel('E');
ylabel('CDF');
title(strcat('CDF for E'));
pause; close

% Forward solve
for x = 1:5000
    y = rand;
    [~,pos] = max(f>y);
    E(x) = xi(pos);
%     plot([0 600],[0 600/E(x)],'Color',blue); hold on
end
% ylabel('strain (mm/mm)')
% xlabel('stress (MPa)')
% title('Forward Solution Samples')
% saveas(gcf,'leforwardsampml','png')
% pause; clf;

maxE = max(E);
for x = 0:10:600
    ind = x/10+1;
    avg(ind) = mean(x./E);
    stdev(ind) = 2*std(x./E);
    minstrain(ind) = min(avg(ind)-x/maxE,stdev(ind));
end
fill([0:10:600 600:-10:0]',[(avg-minstrain)'; flipud((avg+stdev)')],[.25 .25 .25],...
    'Linestyle','none'); hold on
alpha(.5)
plot([0 600],[0 avg(end)],'r','Linewidth',2);
ylabel('strain (mm/mm)')
xlabel('stress (MPa)')
title('Forward Solution Error')

steeldata = importdata('steeldata.txt');% Import data
stress = steeldata(1:201,1);                % Stress is first column
strainavg = mean(steeldata(1:201,2:end),2,'omitnan');
strainstd = std(steeldata(1:201,2:end),0,2,'omitnan'); % Standard deviation of strain
mindiff = min(strainavg-min(steeldata(1:201,2:end),[],2),strainstd);
errorbar(stress, strainavg, mindiff, strainstd,'.','Color',blue);
legend('Bayesian Forward Solve','Bayesian Average','Experimental Data');

saveas(gcf,'leforwarderrml','png')
pause; close 

