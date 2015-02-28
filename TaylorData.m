% Data Matrix
% Impact of meditation training on the default mode network during a restful state
% Matthew Johnson 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Labels
Labels={'PC/PCC', 'DMPFC', 'VMPFC', 'R IPL', 'L IPL', 'R ITC', 'L ITC' ,'R PHG', 'L PHG'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coordinates for seed regions within DMN
Locations= [8 -53 27;
            -10 57 19;
            -2 47 -10;
            48 -56 27;
            -40 -67 34;
            56 -4 -22;
            -56 -10 -17;
            26 -32 -15;
            -26 -29 -18];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Plot brain regions
figure
scatter3(Locations(:,1),Locations(:,2),Locations(:,3));
dx = 3; dy = 3; dz = 3; % displacement so the text does not overlay the data points
text(Locations(:,1)+dx, Locations(:,2)+dy, Locations(:,3)+dz, Labels');
title('Coordinates for seed regions within DMN');
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis') %Label Axes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Correlation Matrix
% Experienced 
EC=[0.00 0.41 0.44 0.53 0.54 0.29 0.29 0.24 0.27;
    0 0.00 0.53 0.25 0.54 0.29 0.29 0.24 0.27;
    0 0 0.00 0.25 0.35 0.39 0.33 0.06 0.11;
    0 0 0 0.00 0.44 0.26 0.23 0.17 0.16;
    0 0 0 0 0.00 0.24 0.26 0.19 0.24;
    0 0 0 0 0 0.00 0.52 0.18 0.19;
    0 0 0 0 0 0 0.00 0.12 0.25;
    0 0 0 0 0 0 0 0.00 0.51;
    0 0 0 0 0 0 0 0 0.00;];
%Inexperienced
IC=[0.00 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.37 0.00 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.45 0.40 0.00 0.0 0.0 0.0 0.0 0.0 0.0;
    0.65 0.34 0.32 0.00 0.0 0.0 0.0 0.0 0.0;
    0.47 0.29 0.29 0.53 0.00 0.0 0.0 0.0 0.0;
    0.31 0.18 0.27 0.27 0.13 0.00 0.0 0.0 0.0;
    0.32 0.29 0.27 0.27 0.17 0.53 0.00 0.0 0.0;
    0.22 -0.10 -0.02 0.16 0.12 0.20 0.14 0.00 0.0;
    0.23 -0.01 0.01 0.17 0.16 0.22 0.27 0.54 0.00;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the Correlation with respect to the Locations 

wgPlot(EC,Locations)
Labels={'PC/PCC', 'DMPFC', 'VMPFC', 'RIPL', 'LIPL', 'RTC', 'LTC' ,'R PHG', 'L PHG'};
text(Locations(:,1)+1, Locations(:,2)+1, Locations(:,3)+1, Labels');
title('Experienced Correlations')
grid on

wgPlot(IC,Locations)
Labels={'PC/PCC', 'DMPFC', 'VMPFC', 'RIPL', 'LIPL', 'RTC', 'LTC' ,'R PHG', 'L PHG'};
text(Locations(:,1)+1, Locations(:,2)+1, Locations(:,3)+1, Labels');
title('Inexperienced Correlations')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Correlation Plots
%Plot Experienced Correlation Matrix
TransposeEC=(EC)';
EC= EC+TransposeEC;
imagesc(EC)
set(gca,'XLim',[0.5 9.5], 'YLim',[0.5 9.5])
set(gca,'XTick',[1:9], 'YTick',[1:9])  
set(gca,'XTickLabel',Labels, 'YTickLabel',Labels)
title('Experienced Correlation Matrix')

%Plot Inexperienced Correlation Matrix
TransposeIC=(IC)';
IC= IC+TransposeIC;
imagesc(IC)
set(gca,'XLim',[0.5 9.5], 'YLim',[0.5 9.5])
set(gca,'XTick',[1:9], 'YTick',[1:9])  
set(gca,'XTickLabel',Labels, 'YTickLabel',Labels)
title('Inxperienced Correlation Matrix')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating axonal volume
Distances=zeros(length(Locations(:,1)),length(Locations(:,1)));
ECVolume=Distances;
ICVolume=Distances;
for i=1:length(Locations(:,1))
    for j=1:length(Locations(:,1))
        Distances(i,j)=sqrt((Locations(i,1)-Locations(j,1))^2+...
            (Locations(i,2)-Locations(j,2))^2)+...
            (Locations(i,3)-Locations(j,3))^2;
        ECVolume(i,j)=Distances(i,j)*EC(i,j);
        ICVolume(i,j)=Distances(i,j)*IC(i,j);
    end
end

TotalOrigECVol=sum(sum(ECVolume));
TotalOrigICVol=sum(sum(ICVolume));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine all possible combinations of brain region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Experienced
Areas=perms(1:9);
Locs=Locations;

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
NewVolumes=NewVolumes./TotalOrigECVol;
plot(sort(NewVolumes))
xlabel('Arrangement Numbers');
ylabel('Amount of Neural Activity compared to original');
title('Experienced Correlation: All Possible Arrangements of Brain Regions');
axis([0 362881 0.6 1.5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine all possible combinations of brain region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inexperienced
Areas=perms(1:9);
Locs=Locations;

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
    z=find(IC~=0);
    [j k]=find(IC~=0);
    NewVolumes(i,1)=sum((sqrt((Locs(j,1)-Locs(k,1)).^2+(Locs(j,2)-Locs(k,2)).^2)+...
        (Locs(j,3)-Locs(k,3)).^2).*IC(z));
end
NewVolumes=NewVolumes./TotalOrigICVol;
plot(sort(NewVolumes))
xlabel('Arrangement Numbers');
ylabel('Amount of Neural Activity compared to original');
title('Inexperienced Correlation: All Possible Arrangements of Brain Regions');
axis([0 362881 0.6 1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partial Correlation Matrix
% Experienced and Inexperienced Partial Correlation Matrix
%  EPM=[0.00 0.13 0.21 0.35 0.27 -0.002 0.03 0.11 0.09;
%      0.11 0.00 0.36 -0.001 0.21 -0.013 0.15 -0.18 -0.01;
%      0.27 0.23 0.00 -0.05 0.04 0.22 0.05 0.001 -0.03;
%      0.47 0.09 -0.05 0.00 0.20 0.10 0.03 0.03 -0.04;
%      0.14 0.09 0.09 0.32 0.00 -0.003 0.02 0.05 0.09;
%      0.05 -0.03 0.13 0.09 -0.07 0.00 0.42 0.10 -0.002;
%      0.06 0.17 0.06 0.03 -0.03 0.43 0.00 -0.05 0.16;
%      0.12 -0.14 -0.06 0.02 0.02 0.09 -0.04 0.00 0.47;
%      0.07 -0.04 -0.08 -0.02 0.07 0.04 0.18 0.48 0.00;];

EP =[0.00 0.13 0.21 0.35 0.27 -0.002 0.03 0.11 0.09;
     0.0 0.00 0.36 -0.001 0.21 -0.013 0.15 -0.18 -0.01;
     0.0 0.0 0.00 -0.05 0.04 0.22 0.05 0.001 -0.03;
     0.0 0.0 0.0 0.00 0.20 0.10 0.03 0.03 -0.04;
     0.0 0.0 0.0 0.0 0.00 -0.003 0.02 0.05 0.09;
     0.0 0.0 0.0 0.0 0.0 0.00 0.42 0.10 -0.002;
     0.0 0.0 0.0 0.0 0.0 0.0 0.00 -0.05 0.16;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00 0.47;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00;];
 
IP =[0.00 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.11 0.00 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.27 0.23 0.00 0.0 0.0 0.0 0.0 0.0 0.0;
     0.47 0.09 -0.05 0.00 0.0 0.0 0.0 0.0 0.0;
     0.14 0.09 0.09 0.32 0.00 0.0 0.0 0.0 0.0;
     0.05 -0.03 0.13 0.09 -0.07 0.00 0.0 0.0 0.0;
     0.06 0.17 0.06 0.03 -0.03 0.43 0.00 0.0 0.0;
     0.12 -0.14 -0.06 0.02 0.02 0.09 -0.04 0.00 0.0;
     0.07 -0.04 -0.08 -0.02 0.07 0.04 0.18 0.48 0.00;];
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Partial Correlation Plots
%Plot Experienced Correlation Matrix
TransposeEP=(EP)';
EP= EP+TransposeEP;
imagesc(EP)
set(gca,'XLim',[0.5 9.5], 'YLim',[0.5 9.5])
set(gca,'XTick',[1:9], 'YTick',[1:9])  
set(gca,'XTickLabel',Labels, 'YTickLabel',Labels)
title('Experienced Partial Correlation Matrix')
colorbar

%Plot Inexperienced Correlation Matrix
TransposeIP=(IP)';
IP= IP+TransposeIP;
imagesc(IP)
set(gca,'XLim',[0.5 9.5], 'YLim',[0.5 9.5])
set(gca,'XTick',[1:9], 'YTick',[1:9])  
set(gca,'XTickLabel',Labels, 'YTickLabel',Labels)
title('Inxperienced Partial Correlation Matrix')
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating axonal volume
Distances=zeros(length(Locations(:,1)),length(Locations(:,1)));
EPVolume=Distances;
IPVolume=Distances;
for i=1:length(Locations(:,1))
    for j=1:length(Locations(:,1))
        Distances(i,j)=sqrt((Locations(i,1)-Locations(j,1))^2+...
            (Locations(i,2)-Locations(j,2))^2)+...
            (Locations(i,3)-Locations(j,3))^2;
        EPVolume(i,j)=Distances(i,j)*EP(i,j);
        IPVolume(i,j)=Distances(i,j)*IP(i,j);
    end
end

TotalOrigEPVol=sum(sum(EPVolume));
TotalOrigIPVol=sum(sum(IPVolume));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine all possible combinations of brain region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Experienced
Areas=perms(1:9);
Locs=Locations;

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
    z=find(EP~=0);
    [j k]=find(EP~=0);
    NewVolumes(i,1)=sum((sqrt((Locs(j,1)-Locs(k,1)).^2+(Locs(j,2)-Locs(k,2)).^2)+...
        (Locs(j,3)-Locs(k,3)).^2).*EP(z));
end
NewVolumes=NewVolumes./TotalOrigEPVol;
plot(sort(NewVolumes))
xlabel('Arrangement Numbers');
ylabel('Amount of Neural Activity compared to original');
title('Experienced Partial Correlation: All Possible Arrangements of Brain Regions');
axis([0 362881 0 5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examine all possible combinations of brain region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inexperienced
Areas=perms(1:9);
Locs=Locations;

NewAreas=zeros(362880,3);
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
    z=find(IP~=0);
    [j k]=find(IP~=0);
    NewVolumes(i,1)=sum((sqrt((Locs(j,1)-Locs(k,1)).^2+(Locs(j,2)-Locs(k,2)).^2)+...
        (Locs(j,3)-Locs(k,3)).^2).*IP(z));
end
NewVolumes=NewVolumes./TotalOrigIPVol;
plot(sort(NewVolumes))
xlabel('Arrangement Numbers');
ylabel('Amount of Neural Activity compared to original');
title('Inexperienced Partial Correlation: All Possible Arrangements of Brain Regions');
axis([0 362881 0 5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Pair-wise Relationships
%Correlation values for all pair-wise relationships (Figure 3)
ECPairwise=nonzeros(EC)
hold on
scatter(1:36,ECPairwise,'v')
hold on; % make sure no new plot window is created on every plot command
ICPairwise=nonzeros(IC)
scatter(1:36,ICPairwise,'s')
plot([0 36], [0 0], 'k-'); % plot the horizontal line
ylabel('Correlation Values')
%xlabel('Pairwise Relationships between default network regions')
legend('Experienced Meditators', 'Inexperienced Meditators')
set(gca, 'XTick', []);
axis ([0 36 -0.2 0.8])

%Partial correlation values for all pair-wise relationships (Figure 5)
EPPairwise=nonzeros(EP)
hold on
scatter(1:36,EPPairwise,'v')
hold on; % make sure no new plot window is created on every plot command
IPPairwise=nonzeros(IP)
scatter(1:36,IPPairwise,'s')
plot([0 36], [0 0], 'k-'); % plot the horizontal line
ylabel('Correlation Values')
%xlabel('Pairwise Relationships between default network regions')
legend('Experienced Meditators', 'Inexperienced Meditators')
set(gca, 'XTick', []);
axis ([0 36 -0.3 0.6])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis

[ECpathlength, ECglobeff]=charpath(EC)
[ICpathlength, ICglobeff]=charpath(IC)
[EPpathlength, EPglobeff]=charpath(EP)
[IPpathlength, IPglobeff]=charpath(IP)