%%%%% In the name of allah
% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article titled: "Frontotemporal Network Contribution to Occluded Face Processing"
% Tracking #: 2024-07457R

%%               load data

clc
clear
close all
load ('data.mat')

%%                 Parameter
[Res_z,cond_z]=all_conditions_zebel;
[Res_j,cond_j]=all_conditions_jenab;

no_fa_nt=find(Res_j(:,1)<6 & Res_j(:,4)>0)';
nt_=reshape(no_fa_nt,8,[]);
on_=311:350;
on_fa_nt=reshape(on_,8,[]);
le_oc=[nt_,on_fa_nt];

lev_occ=[];
for j=[1,4,3,7,2,6,5,8]   %% Order of level of occlusion
    lev_occ=[lev_occ;le_oc(j,:)];
end
le_ix=lev_occ';

le_oc= [reshape(le_ix(:,1:2),[],1),reshape(le_ix(:,3:4),[],1), reshape(le_ix(:,5:6),[],1),reshape(le_ix(:,7:8),[],1) ];

%% Mean responce on condition
ty_ix=[];
ty_ix=le_oc;

it_2m=[];
for ss=1:size(it,1)
    for ind=1:size(ty_ix,2)
        it_2m (ss,:,ind) =nanmean(it{ss}(ismember(cond_j,ty_ix(:,ind)),:));
    end
end
pfc_2m=[];
for ss=1:size(pfc,1)
    for ind=1:size(ty_ix,2)
        pfc_2m (ss,:,ind) =nanmean(pfc{ss}(ismember(cond_j,ty_ix(:,ind)),:));
    end
end
%% Normalize 
bas_li=450:500;

n_it=[];
n_it=        (it_2m-repmat(nanmean(it_2m(:,bas_li,:),2),1,size(it_2m,2)))...
    ./(repmat(max(max(it_2m(:,:,:),[],3),[],2),1,size(it_2m,2),size(it_2m,3))...
    -repmat(min(min(it_2m(:,:,:),[],3),[],2),1,size(it_2m,2),size(it_2m,3)));

n_pfc=[];
n_pfc=        (pfc_2m-repmat(nanmean(pfc_2m(:,bas_li,:),2),1,size(pfc_2m,2)))...
    ./(repmat(max(max(pfc_2m(:,:,:),[],3),[],2),1,size(pfc_2m,2),size(pfc_2m,3))...
    -repmat(min(min(pfc_2m(:,:,:),[],3),[],2),1,size(pfc_2m,2),size(pfc_2m,3)));

%% IT Example Neurons
color_fig=[   [0, 0.4470 ,0.7410];
    [0.8500, 0.3250, 0.0980]
    [0.9290, 0.6940, 0.1250]
    [0.4940, 0.1840, 0.5560]
    [0.4660, 0.6740, 0.1880]
    [0.3010, 0.7450, 0.9330]
    [0.6350, 0.0780, 0.1840]]	;

% best_it=[42,46,37,12,63,130,128,127,125,123,122,116,113,112,111,107,106];
% best_it_1=[12,130,111,37,47,44]; % one peak
% best_it_2=[142,104,80,51,53,42]; % Two peak
it_2=it_2m-mean(it_2m(:,bas_li,:),2);

t_pt=400:1100;
t=t_pt-500;

close all
k=0;
best_it_1=[12]; % Example Neuron 12  & 49
for ii=1:size(best_it_1,2)
    figure('Position', [600 500 300 300])
    for j=1:4
        plot(t,smooth(it_2(best_it_1(ii),t_pt,j)*1000,20),"Color",color_fig(j,:),"LineWidth",1.5)
        hold on
    end
    axis tight
    set(gca,'box','off','TickDir','out','LineWidth', 1.2)
    %     set(gca, 'fontsize', 12, 'fontweight', 'bold');
    set(gca, 'fontsize', 12);

    line([0 0],ylim,'Linewidth',0.3,'Color','k','LineStyle','--')
    line([350 350],ylim,'Linewidth',0.3,'Color','k','LineStyle','--')
    xlabel('Time from Stim. Onset (ms)')
    ylabel('Firing rate (Hz)')
    title(['IT # ',num2str(best_it_1(ii))],FontSize=12)
    xlim([-50 500 ])
    xticks([0 150 350 500])
    %     yticks([0 round(max(ylim)/3) round(max(ylim)/2) round(max(ylim))])
    ylim([min(ylim) max(ylim)+5])
    [~, hobj, ~, ~] = legend({'Intact','Low','High','Full'},'Fontsize',12,'Location','Northeast');
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',2);
    legend boxoff
end

%% PFC Example Neurons
color_fig=[   [0, 0.4470 ,0.7410];
    [0.8500, 0.3250, 0.0980]
    [0.9290, 0.6940, 0.1250]
    [0.4940, 0.1840, 0.5560]
    [0.4660, 0.6740, 0.1880]
    [0.3010, 0.7450, 0.9330]
    [0.6350, 0.0780, 0.1840]]	;

best_it=cat_neuron.n_it_2peak_nice;
% best_it=[42,46,37,12,63,130,128,127,125,123,122,116,113,112,111,107,106];
best_it_1=[12,130,111,37,47,44]; % one peak
best_it_2=[142,104,80,51,53,42]; % Two peak
% best_pfc=[30,26,112,111,153];

pfc_2=pfc_2m-mean(pfc_2m(:,450:500,:),2);
best_pfc=[26,30,42,75,118,112,111,115,118,127,137,142,144,145,149,153,161];

t_pt=400:1100;
t=t_pt-500;
close all
k=0;
best_pfc=[26,149]; % Example Neuron
for ii=1:size(best_pfc,2)
    figure('Position', [600 500 300 300])
    for j=1:4
        plot(t,smooth(pfc_2(best_pfc(ii),t_pt,j)*1000,20),"Color",color_fig(j,:),"LineWidth",2)
        hold on
    end
    axis tight
    set(gca,'box','off','TickDir','out','LineWidth', 1.2)
    %     set(gca, 'fontsize', 12, 'fontweight', 'bold');
    set(gca, 'fontsize', 12);
    line([0 0],ylim,'Linewidth',0.3,'Color','k','LineStyle','--')
    line([350 350],ylim,'Linewidth',0.3,'Color','k','LineStyle','--')
    xlabel('Time from Stim. Onset (ms)')
    ylabel('Firing rate (Hz)')
    title(['PFC # ',num2str(best_pfc(ii))])
    xlim([-50 500 ])
    xticks([0 150 350 500])

    %     yticks([0 round(max(ylim)/3) round(max(ylim)/2) round(max(ylim))])

    ylim([min(ylim) max(ylim)+5])
    [~, hobj, ~, ~] = legend({'Intact','Low','High','Full'},'Fontsize',12,'Location','Northeast');
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',2.5);

    legend boxoff
end







