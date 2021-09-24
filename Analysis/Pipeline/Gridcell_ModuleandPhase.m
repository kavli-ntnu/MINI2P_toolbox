Binsize=12;
SpeedThreadhold=2.5;
%% select all cells with SNR>3, events>100 for analysis

% pick up cells in NeuronMatrix
% colum1: FOV ID;
% colum2: Cell ID;
% colum3: stitched x position
% colum4: stitched y position
% colum5: plane
% colum6: belong to MEC or not
% colum7: SNR
% colum8: event count
% colum9: Is grid cell or not
% colum10: grid score;
% colum11: Is repeated cell or not;
% colum12: Best event train shift;
% colum13: Filted grid cell;
% colum14: Is HD cell or not
% colum15: MVL;
% colum16: Is conjunctive grid cell or not;
SelectedCell=[];
TuningMap_all={};
k=1;
j=1;
for Cell_ID=1:1:size(NeuronMatrix,1)

        FOV_ID=NeuronMatrix(Cell_ID,1);
        i=NeuronMatrix(Cell_ID,2);
        Bestshift=NeuronMatrix(Cell_ID,12);
        SelectedFrame_raw=find(~isnan(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+10)));
        SelectedFrame_filtered=intersect(find(~isnan(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+10))),find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        SelectedFrame_filtered=intersect(SelectedFrame_filtered,find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadhold
        Event_filtered=intersect(find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4*i+12)>0),find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,6)==1));% filter out the frames with speed valid
        Event_filtered=intersect(Event_filtered,find(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,5)>SpeedThreadhold));% filter out the frames with speed threadholds       
        if length(Event_filtered)>100 && StichingPoor{FOV_ID,1}.Information.ExperimentInformation.CellSNR(i)>3 && ~ismember(i,StichingPoor{FOV_ID,1}.Information.ExperimentInformation.RepeatCell) 
            SelectedCell(k)=Cell_ID;
            EventTrain=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(SelectedFrame_filtered,[1 4*i+12]);
            PositionTrain_raw=StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,1:3);
%             PositionTrain_raw(:,2)= PositionTrain_raw(:,2)+(Bestshift* cos(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4) * pi/180));
%             PositionTrain_raw(:,3)= PositionTrain_raw(:,3)+(Bestshift* sin(StichingPoor{FOV_ID,1}.NAT.NAT{1,j}(:,4) * pi/180));
            PositionTrain_try=PositionTrain_raw(SelectedFrame_filtered,:)       ;
            TuningMap=SpatialTuning_BNT.map(PositionTrain_try,EventTrain,'smooth',0,'binWidth',Binsize,'minTime',0.1,'limits',[-48 48 -48 48]);
            TuningMap_all{k,1}=TuningMap.z;
            k=k+1;
        else
        end
Cell_ID
end
%% Calculate autocorrelograms
rmsz = size(TuningMap_all{1,1});
acsz = (rmsz*2 - 1);
nac = prod(acsz);
nu=length(SelectedCell);
acAll = zeros(nac, nu);
for u = 1:nu
    z = TuningMap_all{u,1};
    % interpolate nan values in rate map
    z = inpaintn(z, 10);
    % Use the same autocorrelation method as BNT acorr, but without
    % cropping the size of the output (we'll do this later)
    zac = normxcorr2_general(z, z);
    % Store each autocorrelogram in this matrix (bins x cells)
    acAll(:, u) = zac(:);
end
%%
% Calculate which bins in the autocorrelogram are 'valid' for this
% analysis.
%
% Excluding the center field seems to help, as does excluding the bins
% furthest from the center (i.e. in the corners).

cen = (acsz+1)/2;                          % center position
minr = 3;                                 % minimum radius
maxr = min(rmsz) * 1;                    % maximum radius
[rr,cc] = ndgrid(1:acsz(1), 1:acsz(2));
rr = rr-cen(1);
cc = cc-cen(2);
cendist = hypot(rr, cc);
vrad = cendist < maxr & cendist > minr;
z = nan(acsz);
z(vrad) = 1;
close all
figure
imagesc(z)
axis image off
title("Acorr valid bins");
%% Plot some example autocorrelograms
sply = 4;
splx = 6;
nplt = sply*splx;
figure();
unitInds = randperm(nu, nplt);

for n = 1:nplt
    subplot(sply, splx, n);
    iu = unitInds(n);
    z = reshape(acAll(:, iu), acsz);
    alpha = ones(size(z));
    alpha(~vrad) = 0.5;
    imagesc(z, "alphaData", alpha);
end

axs = findobj(gcf, "type", "axes");
axis(axs, "off", "image", "xy");

%% Do module clustering with UMAP

X = acAll(vrad, :)';
% X = zscore(X);
% Check X has been calculated correctly
% if ~isequal(X, dat.X)
%     warning("Contents of 'X' do not match the pre-calculated values in the sample data files");
% end

% set the random seed to defaut value for reproducibility
rng("default");
warning("off", "MATLAB:DELETE:FileNotFound");
[XUmap, umapObj] = run_umap(X, ...
    "metric", "manhattan", ...
    "n_neighbors",10, ...              % for smaller datasets, you might want to reduce this
    "min_dist", 0.05+eps, ...
    "method", 'c++');
warning("on", "MATLAB:DELETE:FileNotFound");

% Check that UMAP found the same result as I did (it should be
% % reproducible if you use the same implementation and parameters)
% if ~isequal(XUmap, dat.XUmap)
%     warning("Contents of 'XUmap' do not match the pre-calculated values in the sample data files");
% end

% Use DBSCAN to cluster the points in the 2-D UMAP output
mcluIds = dbscan(XUmap, 0.65, 20);
% if ~isequal(mcluIds, dat.mcluIds)
%     warning("Contents of 'mcluIds' do not match the pre-calculated values in the sample data files");
% end

% Check how many clusters we have
mcluIdsU = unique(mcluIds);
nclu = numel(mcluIdsU);
CMP=hsv(nclu+1);
fprintf("Found %u clusters\n", nclu);
 close all
 figure
 scatter(XUmap(:,1),XUmap(:,2),20,CMP(2+mcluIds,:),'filled')
 set(gca,'color',[0 0 0]);
set(gcf,'color',[0 0 0]);
% % Save results of this umap run for comparison
% fn = sprintf("umap_res_%s.mat", datestr(now(), "yyyy-mm-dd_HH-MM-SS"));
% save(fn, "X", "XUmap", "mcluIds", "nclu", "-v7.3", "-nocompression");
%
% Gather information about each cluster
clear moduleClu

for c = 1:nclu

    clear mclu

    % Get indices of all units in the current module cluster
    id = mcluIdsU(c);
    v = mcluIds==id;
    iu = find(v);

    % Store unit identity info in module struct
    mclu.id = id;
    mclu.Cell_indices = iu;
    
    ID=SelectedCell(iu);
    
    mclu.FOV_ID = NeuronMatrix(ID,1);
    mclu.Cell_ID = NeuronMatrix(ID,2);
    mclu.n_total = sum(v);
    mclu.acorrs = reshape(acAll(:, v), acsz(1), acsz(2), []);

    for u = 1:mclu.n_total
        k=mclu.Cell_indices(u);
        mclu.TuningMaps(:, :, u) = TuningMap_all{k,1};
%         mclu.ratemaps_hq(:, :, u) = unit.tc.posFine.z;
%         mclu.tc_hd(:, u) = unit.tc.hd.z;
%         mclu.tc_theta(:, u) = unit.tc.thetaPhase.z;
    end

    % Infer whether what cluster type this is (possibilities are "noise",
    % "non-grid", or "grid"). The DBSCAN output is predictable enough to
    % infer this from the cluster ID.
    %
    % "noise":    Always assigned to ID -1, but not present in every
    %             clustering result.
    %
    % "non-grid": In a good outcome, this should be the single largest
    %             cluster. In every such case which I've seen, this cluster
    %             has been assigned to ID 1.
    %
    % "grid":    In a good outcome, all cluster IDs greater than 1 should
    %            correspond to grid modules.

    % Ascertain the type of cluster, based on the rules/assumptions above
    mclu.name=['cluster ', num2str(mclu.id)];
%     if mclu.id == -1
%         mclu.name = "cluster ";
%     elseif mclu.id == 1
%         mclu.name = "non-grid";
%     elseif mclu.id > 1
%         mclu.name = "grid";
%     end
    moduleClu(c) = mclu;
end



%% plot  whole-session activity may and autocorrelation maps of calcium activity for each cluster
SelectCluster=4
% close all

figure
x0=10;
y0=10;
width=3400;
height=2000;
set(gcf,'position',[x0,y0,width,height])
Column=10;
CelltoCheck=200;
Start=0;
CelltoCheck=moduleClu(1,SelectCluster).n_total;
% for i=1:1:moduleClu(1,SelectCluster).n_total
for i=1:1:CelltoCheck  
        subplot(Column,ceil(CelltoCheck/Column),i,'align')

%     subplot(Column,ceil(moduleClu(1,SelectCluster).n_total/Column),i,'align')
    MAP1=moduleClu(1,SelectCluster).TuningMaps(:,:,i+Start);
    MAP1(isnan(MAP1))=0;
    MAP1=MAP1./prctile(MAP1(:),100);
    MAP2=moduleClu(1,SelectCluster).acorrs(:,:,i+Start);
    MAP2=(MAP2-min(MAP2(:)))./(1-min(MAP2(:)));
    MAP1=imresize(MAP1,size(MAP2),'nearest');
    MAP_combine=[MAP1,MAP2];
    imagesc(flipud(MAP_combine));
    CMP=WJplots.CMP.inferno(256);
    %             colormap(CMP)
    colormap(jet)
    caxis([0 1] );
    ylim([0 size(MAP_combine,1)])
    xlim([0 size(MAP_combine,2)])
%     title(['#',num2str(GridCellAnalysis.IsGridCell{j,1}(i)),' GC:',num2str(GridCellAnalysis.GridScore_shuffled(GridCellAnalysis.IsGridCell{j,1}(i),Shuffling+3,j),'%.2f'),' P:',num2str(max(max(MAP1)),'%.2f')]);
%     title(['#',num2str(GridCellAnalysis.IsGridCell{j,1}(i))]);
    daspect([1 1 1]);
    box off
    axis off
end    
    
    
%%













% calculate mean grid spacing for each grid cluster
% gscores = arrayfun(@(u) u.gridStats.score, U);
% gstats = [U.gridStats]';
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To make things easier, re-order the grid clusters according to their
% grid spacing

% Calculate average grid spacing for each module cluster
meanGspacing = arrayfun(@(s) nanmean(s.spacing), gstats);
isgridclu = [moduleClu.id] >= 2;
moduleCluGrid = moduleClu(isgridclu);
ngridclu = sum(isgridclu);
clear cluGridSpacing
cluGridSpacing = [];
for c = 1:ngridclu
    mclu = moduleCluGrid(c);
    inds = mclu.unit_indices;
    cluGridSpacing(c) = nanmean(meanGspacing(inds));
end

% sort grid clusters by spacing and rename
[~, isort] = sort(cluGridSpacing);
moduleCluGrid = moduleCluGrid(isort);
for c = 1:ngridclu
    moduleCluGrid(c).name = sprintf("grid M%u", c);
    moduleCluGrid(c).id = c+1;
end
moduleClu = [moduleClu(~isgridclu), moduleCluGrid];

mcluIds = nan(nu, 1);
for c = 1:nclu
    mclu = moduleClu(c);
    mcluIds(mclu.unit_indices) = mclu.id;
    fprintf("cluster ID %d, '%s', n=%u units\n", mclu.id, mclu.name, mclu.n_units);
end
disp(' ');




