% Data 
% Meditation-State Functional Connectivity (msFC): Strengthening 
% of the Dorsal Attention Network and Beyond
% Aurthor: Matthew Johnson 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Between Group Effects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default Mode Network

% Labels
DMNLabels={'Posterior cingulate/precuneus', 'Medial Prefrontal', ...
    'Left Lateral Parietal', 'Right Lateral Parietal', 'Left Inferior Temporal',...
    'Right Inferior Temporal', 'Medial Dorsal Thalamus',...
    'Left Posterior Cerebellum','Right Posterior Cerebellum'};

DMNLocations=[...
    8 -53 27;
    2 52 10;
    -46 -62 33;
    44 -62 31;
    -56 -10 -17;
    56 -4 -22;
    0, -12 4;
    -32 -68 -38;
    18 -72 -48;]
    
DMN=[ 0 0 0 0 0 0 0 0 0;
    0.11355318 0 0 0 0 0 0 0 0;
    0.41683059 0.32061067 0 0 0 0 0 0 0;
    -0.03990566 -0.73575101 1.54346167 0 0 0 0 0 0;
    0.93290475 -0.63394433 0.05268126 0.17409072 0 0 0 0 0;
    0.59102585 -0.65045136 -0.60608072 -0.98823958 1.21528256 0 0 0 0;
    -0.27911105 0.09609345 0.40756388 0.1265915 0.28706152 0.91792825 0 0 0;
    -0.42347087 0.96663747 0.80656664 -0.62830591 1.51431815 1.10110768 -0.04309445 0 0;
    -0.5005088 0.29482851 0.90852144 0.30401145 0.0854342 1.53892394 0.65650075 1.81553551 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Plot DMN
TransposeDMN=(DMN)';
DMN= DMN+TransposeDMN;
imagesc(DMN)
set(gca,'XLim',[0.5 9.5], 'YLim',[0.5 9.5])
set(gca,'XTick',[1:9], 'YTick',[1:9])  
set(gca, 'YTickLabel',DMNLabels, 'XTickLabel', [])
title('Default Mode Network')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the Correlation with respect to the Locations 

wgPlot(DMN,DMNLocations)
DMNLabels2={'PC/PCC', 'MPFC', 'L LP', 'R LP', 'L IT', 'R IT', 'MDT','L PC','R PC'};
text(DMNLocations(:,1)+1, DMNLocations(:,2)+1, DMNLocations(:,3)+1, DMNLabels2');
title('Default Mode Network Correlations')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating axonal volume
Distances=zeros(length(DMNLocations(:,1)),length(DMNLocations(:,1)));
DMNVolume=Distances;
for i=1:length(DMNLocations(:,1))
    for j=1:length(DMNLocations(:,1))
        Distances(i,j)=sqrt((DMNLocations(i,1)-DMNLocations(j,1))^2+...
            (DMNLocations(i,2)-DMNLocations(j,2))^2)+...
            (DMNLocations(i,3)-DMNLocations(j,3))^2;
        DMNVolume(i,j)=Distances(i,j)*DMN(i,j);
    end
end
TotalDMNVol=sum(sum(DMNVolume));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine all possible combinations of brain region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Experienced
Areas=perms(1:9);
Locs=DMNLocations;

NewVolumes=zeros(362880,1);

for i=1:362880
    Locs(1,:)=Locations(Areas(i,1),:);
    Locs(2,:)=Locations(Areas(i,2),:);
    Locs(3,:)=Locations(Areas(i,3),:);
    Locs(4,:)=Locations(Areas(i,4),:);
    Locs(5,:)=Locations(Areas(i,5),:);
    Locs(6,:)=Locations(Areas(i,6),:);   
    Locs(7,:)=Locations(Areas(i,7),:);    
    Locs(8,:)=Locations(Areas(i,8),:);    
    Locs(9,:)=Locations(Areas(i,9),:);
    z=find(EC~=0);
    [j k]=find(EC~=0);
    NewVolumes(i,1)=sum((sqrt((Locs(j,1)-Locs(k,1)).^2+(Locs(j,2)-Locs(k,2)).^2)+...
        (Locs(j,3)-Locs(k,3)).^2).*EC(z));
end
NewVolumes=NewVolumes./TotalDMNVol;
plot(sort(NewVolumes))
xlabel('Arrangement Numbers');
ylabel('Amount of Neural Activity compared to original');
title('DMN Correlations: All Possible Arrangements of Brain Regions');
axis([0 362881 0.6 1.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dorsal Attention Network

% Labels
DANLabels={'Left Frontal Eyefield', 'Right Frontal Eyefield', ...
    'Left Posterior IPS', 'Right Posterior IPS', 'Left Anterior IPS', ...
    'Right Anterior IPS', 'Left MT', 'Right MT'};
    

DAN=[ 0 0 0 0 0 0 0 0;
    0.31040452 0 0 0 0 0 0 0;
    -0.06895648 -0.91117578 0 0 0 0 0 0; 
    1.48690868 1.23893297 -1.35158566 0 0 0 0 0;
    0.86765425 -0.44317913 0.67430468 0.45047156 0 0 0 0;
    2.216917 1.13267216 -0.41423442 -0.31084837 0.82881631 0 0 0;
    0.88826522 0.73187341 0.8121265 0.33528222 0.96949473 1.15873859 0 0;
    2.40379655 1.92647685 1.45967123 2.6062609 2.89577123 1.28535601 2.40373963 0];
    

TransposeDAN=(DAN)';
DAN= DAN+TransposeDAN;
imagesc(DAN)
set(gca,'XLim',[0.5 8.5], 'YLim',[0.5 8.5])
set(gca,'XTick',[1:8], 'YTick',[1:8])  
set(gca, 'YTickLabel',DANLabels, 'XTickLabel', [])
title('Dorsal Attention Network')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Executive Control Network

% Labels
ECNLabels={'Dorsal Medial PFC', 'Left Anterior PFC', ...
    'Right Anterior PFC', 'Left Superior Parietal', 'Right Superior Parietal'};
    

ECN=[ 0 0 0 0 0;
    0.16841477 0 0 0 0;
    -0.31091113 -0.69439452 0 0 0;
    -0.73965638 -1.03300666 0.449007621 0 0;
    0.52108483 0.55402123 -0.12536615 0.95449354 0];

TransposeECN=ECN';
ECN= ECN+TransposeECN;
imagesc(ECN)
set(gca,'XLim',[0.5 5.5], 'YLim',[0.5 5.5])
set(gca,'XTick',[1:5], 'YTick',[1:5])  
set(gca, 'YTickLabel',ECNLabels, 'XTickLabel', [])
title('Executive Control Network')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saliance Network

% Labels
SNLabels={'Dorsal Anterior Cingulate', 'Left Anterior PFC', ...
    'Right Anterior PFC', 'Left Insula', 'Right Insula', ...
    'Left Lateral Parietal', 'Right Lateral Parietal'};
    

SN=[ 0 0 0 0 0 0 0;
    0.90898511 0 0 0 0 0 0;
    -0.25383124 -0.38113833 0 0 0 0 0;
    0.02684801 -0.75813955 -1.31835793 0 0 0 0;
    0.41976991 -0.84710253 -1.42870933 0.13888497 0 0 0;
    -0.29589018 -1.78208187 -0.27952114 -0.73573779 -0.71776875 0 0;
    -1.11089951 -1.6437564 -1.10015352 -0.36434026 -1.20335322 1.06406012 0];

TransposeSN=SN';
SN= SN+TransposeSN;
imagesc(SN)
set(gca,'XLim',[0.5 7.5], 'YLim',[0.5 7.5])
set(gca,'XTick',[1:7], 'YTick',[1:7])  
set(gca, 'YTickLabel',SNLabels, 'XTickLabel', [])
title('Saliance Network')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis
%Charecteristic Path Length and global efficiency
[DMNcharpathlength, DMNglobaleff] =charpath(DMN);
[DANcharpathlength, DANglobaleff]=charpath(DAN);
[ECNcharpathlength, ECNglobaleff]=charpath(ECN);
[SNcharpathlength, SNglobaleff]=charpath(SN);




