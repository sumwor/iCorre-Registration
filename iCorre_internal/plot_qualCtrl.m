function plot_qualCtrl(path_QC, option_label)

%% plot the quality control figures
%% supported by NoRMCorre

% load the data 
cd(path_QC);
if strcmp(option_label, 'RMC')
    matFile = dir('*_RM_*.mat');
elseif strcmp(option_label,'NRMC')
    matFile = dir('*_NRM_*.mat');
end

for ii = 1:numel(matFile)
    if ii == 1
        data = load(matFile(ii).name);
    else
        % probably not useful, come in handy when the QC varibles were
        % saved in chunks of trials (to ease the memory issue)
        dataNames = fieldnames(data);
        tempData = load(matFile(ii).name);
        for jj = 1:numel(dataNames)
            if ~strcmp(dataNames{jj}, 'options_label')
                data.(dataNames{jj}) = cat(ndims(data.(dataNames{jj})),data.(dataNames{jj}), tempData.(dataNames{jj}));
            end
        end
    end      
end

%% compute and plot metrics
mY = mean(data.mYc, 3); % average frame of raw data (average of averages of each trial)
nnY = min(data.nnYc);
mmY = max(data.mmYc);
T = sum(data.Tc); % total number of frames
cY = cell2mat(data.cYc');
figure;
ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')

if strcmp(data.options_label, 'RMC')  % rigid motion correction
    mM1 = mean(data.mMRc, 3);
    cM1 = cell2mat(data.cMRc');
    ax2 = subplot(2,2,2);imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,2,3); plot(1:T,cY,1:T,cM1); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold') 
    subplot(2,2,4); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
    xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    print(gcf, '-r0', fullfile(path_QC , 'RMC_metrics'), '-dpng'); %png format

elseif strcmp(data.options_label, 'NRMC')
    mM2 = mean(data.mMNRc, 3);
    cM2 = cell2mat(data.cMNRc');
    ax2 = subplot(2,2,2);imagesc(mM2,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,2,3); plot(1:T,cY,1:T,cM2); legend('raw data','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold') 
    subplot(2,2,4); scatter(cY,cM2); hold on; plot([0.9*min(cY),1.05*max(cM2)],[0.9*min(cY),1.05*max(cM2)],'--r'); axis square;
    xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    print(gcf, '-r0', fullfile(path_QC , 'NRMC_metrics'), '-dpng'); %png format

end

%% plot shifts        
if strcmp(data.options_label, 'RMC')  % rigid motion correction
  shifts_r = cell2mat(data.shifts_rc');
  figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')
    print(gcf, '-r0', fullfile(path_QC , 'RMC_shifts'), '-dpng'); %png format
elseif strcmp(data.options_label, 'NRMC')
  %shifts_nr = cell2mat(data.shifts_nrc);
  shifts_x = cell2mat(data.shifts_xc');
  shifts_y = cell2mat(data.shifts_yc');
  figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM2); legend('raw data','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x);  title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y);  title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')
    print(gcf, '-r0', fullfile(path_QC , 'NRMC_shifts'), '-dpng'); %png format
end



