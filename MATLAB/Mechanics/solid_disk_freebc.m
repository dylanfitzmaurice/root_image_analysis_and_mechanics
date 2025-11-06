% solid_disk_freebc.m  — solid circle, traction-free rim + applys force to outer edge
% Minimal pins remove rigid-body modes; stresses reflect a free boundary.
% solid_disk_freebc_v2.m — solid circle, traction-free rim + uniform pressure
% Solid disk: you're applying a traction on a single outer boundary
% All of the pressure work is transmitted through a continuous solid.

clear; close all; clc

%% ---------------- User knobs ----------------
R            = 100;        % radius (model units)
Hmax         = 3;          % target mesh size
plane        = 'strain';   % 'stress' or 'strain'
Eyoung       = 1e3;        % Young's modulus
nu           = 0.3;        % Poisson's ratio
thickness    = 1.0;        % only used for plane stress
p_rim        = -1;         % uniform pressure on the outer rim, acts outward (tensile) — i.e., it inflates/expands the rim.
use_bodyforce= false;      % set true to test body force instead of pressure
bvec         = [0; 0];     % body force per area (if use_bodyforce=true)


%% ------------- 1) Model & circle geometry of radius R -------------
structuralModel = createpde('structural', ['static-plane' plane]);

% decsg geometry for a single circle of radius R centered at (0,0)
% gd: [1; xc; yc; r] for a circle
gd = [1; 0; 0; R];
ns = char('C1'); ns = ns';
sf = 'C1';
dl = decsg(gd, sf, ns);
geometryFromEdges(structuralModel, dl);

%% ------------- 2) Material -------------
if strcmpi(plane,'stress')
    structuralProperties(structuralModel, ...
        'YoungsModulus',Eyoung, 'PoissonsRatio',nu, 'Thickness',thickness);
else
    structuralProperties(structuralModel, ...
        'YoungsModulus',Eyoung, 'PoissonsRatio',nu);
end

%% ------------- 3) Loads (choose ONE) -------------
if ~use_bodyforce
    % (A) Uniform pressure on the outer boundary (self-equilibrated).
    % Positive Pressure acts inward along the outward normal.
    structuralBoundaryLoad(structuralModel, ...
        'Edge', 1:structuralModel.Geometry.NumEdges, ...
        'Pressure', p_rim);
else
    % (B) Uniform body force (needs pins; not self-equilibrated).
    structuralBodyLoad(structuralModel, 'GravitationalAcceleration', bvec);
end

%% ------------- 4) Mesh -------------
generateMesh(structuralModel, 'Hmax', Hmax);

%% ------------- 5) Minimal constraints (remove rigid-body modes only) -------------
% Pick two geometry vertices on the rim
V = structuralModel.Geometry.Vertices;   % 2 x Nv array of vertex coordinates

[~, vA] = max(V(1,:));         % rightmost vertex
[~, vB] = max(V(2,:));         % topmost vertex
if vB == vA
    [~, vB] = min(V(1,:));     % fall back to leftmost if they coincide
end

% Pin one vertex completely (kills translations) …
structuralBC(structuralModel, 'Vertex', vA, ...
    'XDisplacement', 0, 'YDisplacement', 0);

% … and pin only one DOF at a second vertex (kills rigid rotation)
structuralBC(structuralModel, 'Vertex', vB, ...
    'YDisplacement', 0);


%% ------------- 6) Solve -------------
Rsol = solve(structuralModel);

% --- after you call Rsol = solve(structuralModel); ---

% Mesh node coordinates (2×N)
nodes = structuralModel.Mesh.Nodes;
x = nodes(1,:)'; 
y = nodes(2,:)';

% Radial / tangential unit vectors at nodes
rr = hypot(x,y);
er = [x./max(rr,eps), y./max(rr,eps)];   % radial unit
i0 = rr < 1e-14;                         % guard the center
er(i0,:) = repmat([1,0], nnz(i0), 1);    % arbitrary radial at r≈0
et = [-er(:,2), er(:,1)];                % tangential unit

% Cartesian nodal stresses from solution (N×1 each)
S   = Rsol.Stress;
sxx = S.sxx;
syy = S.syy;
sxy = S.sxy;  % τ_xy

%% ------------- 7) Plots (VM, radial, hoop) -------------
vm = Rsol.VonMisesStress; vmax = max(vm); if vmax<=0, vmax=1; end
fprintf('Solid: mean(VM)/|p| = %.3g\n', mean(vm)/abs(p_rim));

% figure('Color','w');
% pdeplot(structuralModel, 'XYData', vm, ...
%     'Deformation', Rsol.Displacement, 'DeformationScaleFactor', 1, ...
%     'ColorMap', 'parula');
% title(sprintf('Von Mises stress (plane %s), p_{rim}=%.3g', plane, p_rim));
% axis equal tight; caxis([0 vmax]); colorbar

% Cartesian stresses at nodes
S   = Rsol.Stress;  % fields: sxx, syy, sxy (nodal)
x   = nodes(1,:)'; y = nodes(2,:)'; rr = hypot(x,y);
er  = [x./max(rr,eps), y./max(rr,eps)];   % radial unit
et  = [-er(:,2), er(:,1)];                % hoop unit

Sxx = S.sxx; Syy = S.syy; Sxy = S.sxy;
s_rr = er(:,1).*(Sxx.*er(:,1) + Sxy.*er(:,2)) + ...
       er(:,2).*(Sxy.*er(:,1) + Syy.*er(:,2));
s_tt = et(:,1).*(Sxx.*et(:,1) + Sxy.*et(:,2)) + ...
       et(:,2).*(Sxy.*et(:,1) + Syy.*et(:,2));

% Mesh, fields, and a deformed overlay
P   = structuralModel.Mesh.Nodes.';                 % N×2
T3  = structuralModel.Mesh.Elements(1:3,:).';       % Ne×3
ux  = Rsol.Displacement.ux; 
uy  = Rsol.Displacement.uy;
scale = 1 / max(1e-9, sqrt(mean(ux.^2 + uy.^2)));
Pdef  = P + scale * [ux uy];

% If you formed per-element fields, you can pass them too — the helper handles both.
Srr_node = s_rr;                % nodal already (your fprintf checked this)
Stt_node = s_tt;                % nodal
VM_node  = Rsol.VonMisesStress; % nodal

%% If you don't want white at 0, use 'parula' or shift the 'Center'
%plot_tri_scalar(P, T3, Srr_node, 'Title','SOLID — \sigma_{rr}', ...
%                'Pplot',Pdef, 'ShowMesh',true, 'Colormap','parula');

% Radial stress (RWB centered at 0)
Sr = s_rr(:);
S  = max(abs(Sr(isfinite(Sr))));
plot_tri_scalar(P, T3, Sr, 'Pplot', Pdef, 'ShowMesh', true, ...
    'Colormap', 'diverging', 'Center', 0, 'CLim', [-S S]);
title('SOLID — radial stress \sigma_{rr}', 'Interpreter','tex', 'FontWeight','bold');
%\sigma_{rr} / |p|'

% Hoop stress (RWB centered at 0)
St = s_tt(:);
S  = max(abs(St(isfinite(St))));
plot_tri_scalar(P, T3, St, 'Pplot', Pdef, 'ShowMesh', true, ...
    'Colormap', 'diverging', 'Center', 0, 'CLim', [-S S]);
title('SOLID — circumferential stress \sigma_{\theta\theta}', 'Interpreter','tex', 'FontWeight','bold');
%\sigma_{\theta\theta} / |p|'

%% Von Mises is ≥0 → keep sequential (parula)
%plot_tri_scalar(P, T3, Rsol.VonMisesStress, 'Pplot', Pdef, 'ShowMesh', true, ...
%    'Colormap', 'parula');   % no Center/CLim symmetry here

vm = Rsol.VonMisesStress(:);
Vmax = prctile(vm(isfinite(vm)), 99);        % robust cap; or use max(...)
plot_tri_scalar(P, T3, vm, 'Pplot', Pdef, 'ShowMesh', true);
colormap( white_to_red(256) );               % <- defined below
caxis([0 Vmax]); colorbar
title('SOLID — Von Mises stress', 'Interpreter','tex', 'FontWeight','bold');


%% helpers
function triVals = nodalToElement(mesh, nodalVals)
% simple average of nodal values to triangle centroids (for pdeplot XYData)
T = mesh.Elements.'; 
triVals = mean(nodalVals(T), 2);
end

function plot_tri_scalar_bwr(P, T, vals, ttl, clim)
    % Symmetric blue↔white↔red around 0, works with per-node OR per-element vals.
    if nargin < 5 || isempty(clim)
        vmax = max(abs(vals(:))); if vmax<=0, vmax = 1; end
        clim = [-vmax, +vmax];
    end
    [vals_node, Z] = ensure_node_vals(P, T, vals);
    figure('Color','k'); hold on; axis equal off
    trisurf(T, P(:,1), P(:,2), Z, vals_node, 'EdgeColor','k');
    view(2); colormap(bwr_colormap(256)); caxis(clim); colorbar
    title(ttl); drawnow
end

function plot_tri_scalar_bwrpos(P, T, vals, ttl, vmax_in)
    % White at 0 → red for positive (for von-Mises, etc.). Per-node OR per-element.
    if nargin < 5 || isempty(vmax_in)
        vmax_in = max(vals(:)); if vmax_in<=0, vmax_in = 1; end
    end
    [vals_node, Z] = ensure_node_vals(P, T, vals);
    figure('Color','k'); hold on; axis equal off
    trisurf(T, P(:,1), P(:,2), Z, vals_node, 'EdgeColor','none');
    view(2); colormap(bwr_colormap_centered(256, 0, vmax_in, 0)); caxis([0 vmax_in]); colorbar
    title(ttl); drawnow
end

function [vals_node, Z] = ensure_node_vals(P, T, vals)
    % Accept per-node (N) or per-element (Ne) values and return per-node.
    % T may be Ne×nvpe or nvpe×Ne (nvpe = 3 or 6 for triangles).
    nN = size(P,1);

    % Ensure T is Ne×nvpe
    if size(T,2) > size(T,1) && size(T,1) ~= nN
        % likely nvpe×Ne → transpose
        T = T.';
    end
    T = double(T);
    Ne = size(T,1);
    nvpe = size(T,2);

    if numel(vals) == nN
        % Already nodal
        vals_node = vals(:);
    elseif numel(vals) == Ne
        % Element → node: simple equal-weight scatter then average
        idx = T(:);                     % length = Ne*nvpe
        w   = repelem(vals(:), nvpe);   % same length as idx
        sums = accumarray(idx, w, [nN,1], @sum, 0);
        cnts = accumarray(idx, 1, [nN,1], @sum, 0);
        vals_node = sums ./ max(cnts,1);
    else
        error('vals must be length N (nodes) or Ne (elements).');
    end

    Z = zeros(nN,1);  % flat 2D surface
end

function cmap = bwr_colormap_centered(n, vmin, vmax, center)
    if nargin < 1, n = 256; end
    if ~(vmin < center && center < vmax)
        cmap = bwr_colormap(n); return
    end
    f  = (center - vmin) / (vmax - vmin);  % white position in [0,1]
    nL = max(1, round(n * f));
    nR = max(1, n - nL);
    cL = [linspace(0,1,nL)', linspace(0,1,nL)', ones(nL,1)];   % blue→white
    cR = [ones(nR,1), linspace(1,0,nR)', linspace(1,0,nR)'];   % white→red
    cL(end,:) = [1 1 1]; cR(1,:) = [1 1 1];
    cmap = [cL; cR];
end

function [P2,T2,used] = prune_degenerate_with_used(P, T, epsA)
    % Make sure triangles are 3-node connectivity in rows
    if size(T,2) > size(T,1), T = T.'; end
    T = T(:,1:3);

    % Triangle areas
    A = 0.5*abs( (P(T(:,2),1)-P(T(:,1),1)).*(P(T(:,3),2)-P(T(:,1),2)) ...
               - (P(T(:,2),2)-P(T(:,1),2)).*(P(T(:,3),1)-P(T(:,1),1)) );

    if nargin<3 || isempty(epsA), epsA = 1e-12; end
    Ath  = epsA * max(1, median(A(A>0)));
    keep = A > Ath;

    T2   = T(keep,:);
    used = unique(T2(:));
    map  = zeros(size(P,1),1); map(used) = 1:numel(used);
    P2   = P(used,:);
    T2   = reshape(map(T2), size(T2));
end

function vals_node = to_nodal(P, T, vals)
    % Accept per-node (N) or per-element (Ne) and return per-node on (P,T)
    N  = size(P,1);
    if size(T,2) > size(T,1), T = T.'; end
    T  = T(:,1:3);
    Ne = size(T,1);
    if numel(vals) == N
        vals_node = vals(:);
    elseif numel(vals) == Ne
        % average element values to incident nodes
        idx  = T(:);
        w    = repelem(vals(:), 3);
        sums = accumarray(idx, w, [N,1], @sum, 0);
        cnts = accumarray(idx, 1, [N,1], @sum, 0);
        vals_node = sums ./ max(cnts,1);
    else
        error('vals must be length N (nodal) or Ne (elemental).');
    end
end

function cmap = bwr_colormap(n)
    if nargin<1, n=256; end
    n1=floor(n/2); n2=n-n1;
    c1=[linspace(0,1,n1)', linspace(0,1,n1)', ones(n1,1)];
    c2=[ones(n2,1), linspace(1,0,n2)', linspace(1,0,n2)'];
    c1(end,:)=1; c2(1,:)=1;
    cmap=[c1; c2];
end

function cmap = bwr_centered(n, vmin, vmax, center)
    if nargin<1, n=256; end
    if ~(vmin < center && center < vmax), cmap = bwr_colormap(n); return; end
    f  = (center - vmin) / (vmax - vmin);
    nL = max(1, round(n * f)); nR = max(1, n - nL);
    cL = [linspace(0,1,nL)', linspace(0,1,nL)', ones(nL,1)];
    cR = [ones(nR,1), linspace(1,0,nR)', linspace(1,0,nR)'];
    cL(end,:) = [1 1 1]; cR(1,:) = [1 1 1];
    cmap = [cL; cR];
end

function plot_tri_scalar(P, T, vals, varargin)
% plot_tri_scalar(P, T, vals, 'Title',ttl, 'Pplot',Pplot, 'ShowMesh',true, ...
%                 'Colormap','diverging'|'parula', 'Center',0, 'CLim',[])
    opt = struct('Title','', 'Pplot',[], 'ShowMesh',true, ...
                 'Colormap','diverging', 'Center', 0, 'CLim',[]);
    for k=1:2:numel(varargin), opt.(varargin{k}) = varargin{k+1}; end

    % Prune once using ORIGINAL coordinates (keeps indexing consistent)
    [P2, T2, used] = prune_degenerate_with_used(P, T, 1e-12);

    % Values → nodal on the pruned mesh
    if numel(vals) == size(P,1)
        vals_node_full = vals(:);
        vals2 = vals_node_full(used);
    else
        % If elemental, prune the element list first then average to nodes
        if size(T,2) > size(T,1), Tall = T.'; else, Tall = T; end
        Tall = Tall(:,1:3);
        A = 0.5*abs( (P(Tall(:,2),1)-P(Tall(:,1),1)).*(P(Tall(:,3),2)-P(Tall(:,1),2)) ...
                   - (P(Tall(:,2),2)-P(Tall(:,1),2)).*(P(Tall(:,3),1)-P(Tall(:,1),1)) );
        Ath  = 1e-12 * max(1, median(A(A>0)));
        keep = A > Ath;
        vals_keep = vals(keep);
        vals2 = to_nodal(P2, T2, vals_keep);
    end

    % Choose plotting coordinates (original or deformed), using SAME nodes
    if isempty(opt.Pplot), Pplot2 = P2; else, Pplot2 = opt.Pplot(used,:); end

    % Draw
    figure('Color','k'); hold on
    trisurf(T2, Pplot2(:,1), Pplot2(:,2), zeros(size(Pplot2,1),1), vals2, ...
            'EdgeColor','none', 'FaceColor','interp');
    view(2); axis equal tight; colorbar
    if opt.ShowMesh
        triplot(T2, Pplot2(:,1), Pplot2(:,2), 'k-', 'LineWidth', 0.3);
    end
    title(opt.Title);

    % Colormap
    switch lower(opt.Colormap)
        case 'parula'
            colormap(parula);
        case 'diverging'
            % White at "Center" (0 by default). Move Center if you want to avoid white.
            vmin = min(vals2(isfinite(vals2))); vmax = max(vals2(isfinite(vals2)));
            if isempty(opt.CLim)
                caxis([vmin vmax]);
            else
                caxis(opt.CLim);
                vmin = opt.CLim(1); vmax = opt.CLim(2);
            end
            colormap(bwr_centered(256, vmin, vmax, opt.Center));
        otherwise
            colormap(parula);
    end
end

function cmap = white_to_red(n)
    if nargin<1, n=256; end
    r = ones(n,1); g = linspace(1,0,n)'; b = linspace(1,0,n)';  % white→red
    cmap = [r g b];
end
