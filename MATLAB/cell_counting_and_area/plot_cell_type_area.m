%plot_cell_type_area.m

clear; close all; clc

% ---- List your sample files ----
number_of_samples = 6;
files = arrayfun(@(k) sprintf('WT_16DAG_root%d.mat', k), 1:number_of_samples, 'uni', 0);
S = numel(files);

% ---- Error bar mode: 'sem' or 'sd' ----
err_mode = 'sem';
ylab = 'Cell area (\mum^2)';
default_type_names = {'Epidermis','Cortex','Middle Cortex','Endodermis','Stele'};
type_names = [];
T = [];

% ---- Colors ----
cols_type = [0.10 0.45 0.85; 0.10 0.65 0.10; 0.55 0.10 0.70; 0.70 0.40 0.05; 0.50 0.50 0.50];
cols_samp = lines(S);                     % distinct sample colors for dots

% ---- Storage ----
mu_samp   = [];                    % [T x S] mean area per type per sample
n_samp_t  = [];                    % [T x S] counts per type/sample
sample_cells = {};                 % {S,T} raw per-cell vectors by sample & type

% ================== Load & normalize each sample ==================
for si = 1:S
    Sx = load(files{si});

    % Determine type_names / T from first file
    if isempty(type_names)
        if isfield(Sx,'type_names') && ~isempty(Sx.type_names)
            type_names = Sx.type_names(:).';
        else
            type_names = default_type_names;
        end
        T = numel(type_names);
        if size(cols_type,1) ~= T
            cols_type = repmat(cols_type, ceil(T/size(cols_type,1)), 1);
            cols_type = cols_type(1:T,:);
        end
        mu_samp  = NaN(T,S);
        n_samp_t = zeros(T,S);
        sample_cells = cell(S,T);
    end

    % Accept either (A0, cell_type) OR areas{t} from your single-sample saves
    if isfield(Sx,'A0') && isfield(Sx,'cell_type')
        for t = 1:T
            idx = (Sx.cell_type==t) & isfinite(Sx.A0) & Sx.A0>0;
            vals = Sx.A0(idx);
            sample_cells{si,t} = vals(:);
            n_samp_t(t,si) = numel(vals);
            if ~isempty(vals), mu_samp(t,si) = mean(vals); end
        end
    elseif isfield(Sx,'areas') && iscell(Sx.areas)
        for t = 1:min(T,numel(Sx.areas))
            vals = Sx.areas{t};
            vals = vals(isfinite(vals) & vals>0);
            sample_cells{si,t} = vals(:);
            n_samp_t(t,si) = numel(vals);
            if ~isempty(vals), mu_samp(t,si) = mean(vals); end
        end
    else
        error('File %s lacks required variables. Expect A0+cell_type or areas{t}.', files{si});
    end
end

% ================== Aggregate for bars ==================
mu_plot  = mean(mu_samp, 2, 'omitnan');         % bar height = mean of per-sample means
sd_plot  = std(mu_samp, 0, 2, 'omitnan');
n_sampOK = sum(~isnan(mu_samp), 2);             % number of contributing samples per type
sem_plot = sd_plot ./ sqrt(max(1,n_sampOK));
err_plot = strcmpi(err_mode,'sd') .* sd_plot + ~strcmpi(err_mode,'sd') .* sem_plot;

% total cell counts per type across all samples
N_cells_total = sum(n_samp_t, 2);

% ================== Plot ==================
figure('Color','k'); hold on
b = bar(mu_plot, 'FaceColor','flat'); b.CData = cols_type;
errorbar(1:T, mu_plot, err_plot, 'w', 'LineStyle','none', 'LineWidth', 1.2);

% ---- Overlay EVERY cell measurement, colored by sample, offset per sample ----
off = linspace(-0.28, +0.28, S);   % symmetric offsets around each bar
jitter = 0.04;                     % small jitter within each sample's lane
ptSize = 16;

for si = 1:S
    for t = 1:T
        y = sample_cells{si,t};
        if isempty(y), continue; end
        xj = (t + off(si)) + jitter*randn(size(y));
        scatter(xj, y, ptSize, 'filled', ...
            'MarkerFaceColor', cols_samp(si,:), ...
            'MarkerFaceAlpha', 0.35, ...
            'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);
    end
end

% Optional: overlay per-sample means as larger markers (uncomment if wanted)
% for si = 1:S
%     for t = 1:T
%         if ~isnan(mu_samp(t,si))
%             scatter(t + off(si), mu_samp(t,si), 48, ...
%                 'MarkerFaceColor', cols_samp(si,:), 'MarkerEdgeColor','w', 'LineWidth', 0.6);
%         end
%     end
% end

% Annotate counts above bars
for t = 1:T
    ytxt = mu_plot(t) + (err_plot(t) + 1e-12)*1.08;
    text(t, ytxt, sprintf('S=%d | N=%d', n_sampOK(t), N_cells_total(t)), ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
         'FontSize', 9, 'Color','w');
end

% Cosmetics
ax = gca; ax.Color='k'; ax.XColor='w'; ax.YColor='w';
ax.GridColor=[0.6 0.6 0.6]; ax.GridAlpha=0.15;
xlim([0.5, T+0.5]); ylim([0, inf]);
set(ax, 'XTick', 1:T, 'XTickLabel', type_names, 'XTickLabelRotation', 20);
ylabel(ylab, 'Color','w');
title(sprintf('Cell area by cell-type across %d samples (bars = mean of per-sample means, error = %s)', ...
      S, upper(err_mode)), 'Color','w');
grid on; box on

% Legend for samples
lgdLabels = arrayfun(@(k) sprintf('Sample %d',k), 1:S, 'uni', 0);
plot(NaN,NaN,'w-'); % spacer so legend has a first bar entry if desired
%legend(lgdLabels, 'TextColor','w', 'Location','bestoutside', 'Color',[0.1 0.1 0.1]);

% ================== Console summary ==================
fprintf('\nPer-type stats across samples (%s of per-sample means):\n', upper(err_mode));
for t = 1:T
    fprintf('  %-15s  S=%d  N=%d  mean=%.4g  SD=%.4g  SEM=%.4g\n', ...
        type_names{t}, n_sampOK(t), N_cells_total(t), mu_plot(t), sd_plot(t), sem_plot(t));
end
