% cell_type_counting_and_area.m
%
% Generally, this code loads a transverse optical section .svg, the user
% made using Abode Illustrator, then it meshes the cell walls and applies
% a 2x radially thicker outer cell wall and plots this visually. This 
% visualation is useful as if the .svg file needs editing it will show
% internal cells which need to be editted. Next, the code allows the user
% to select cell type by clicking each cell. After cell type selcetion, the
% number of each cell type and area for each cell is calculate.
%
% 
% Workflow:
% User inputs & units
% SVG → closed loops (cells) from <polygon> + <path>
% Calibrate to microns (viewBox shift + per-axis scale)
% Disjoint polyshapes (overlap resolution) → Pcells
% Bands + meshing (shared nodes, constrained DT, inner/exterior edge classification)
% (Optional) Out-of-plane thickness maps for FE weighting
% Radial (in-plane) outer thickening + re-mesh (geometry encodes thickness)
% Sanity viz: added radial strips / perimeter check
% Interactive cell-type assignment (click lumens)
% Per-type area stats & deviation from type mean
% Helper functions (densify, snap map, edges)


clear; close all; clc
%% ====================== 1) USER INPUTS & UNITS ============================
% Purpose:
%   - Point to the SVG to analyze.
%   - Define pixel→micron scale (can be anisotropic if SVG units differ in X/Y).
% Why:
%   - All geometric ops and downstream measurements (areas, thickness) should be in µm.
% Knobs:
%   - PX2UM_X, PX2UM_Y (µm per SVG unit). If your SVG is already in µm, set both = 1.

% [Pick your SVG]
svgFile =  '/Users/dylan/Downloads/SUM_Reslice of Reslice_small_pep_example.svg';

% Set your constant pixel-to-micron scale (square pixels assumed).
PX2UM_X = 0.207569;    % µm per user-unit in X
PX2UM_Y = 0.207569;    % µm per user-unit in Y  (use a different value if anisotropic)

name_for_saving = 'small_peptide_root5';


%% --------------- 2) SVG PARSING → CLOSED LOOPS (candidate cells) ----------
% Purpose:
%   - Read <polygon> and <path> data and convert to explicit (x,y) loops.
%   - For Beziers (C/Q), sample curves into polylines (NSEG points/segment).
% Why:
%   - Downstream MATLAB polyshape/mesh tools need explicit vertices.
% Knobs:
%   - NSEG (curve sampling density); LINE_WIDTH only affects preview plots.
% Pitfalls:
%   - Empty 'points' or 'd' attributes are skipped.
%   - Ensure subpaths are closed (Z) so first==last.
%   - Beware scientific notation in numeric tokens (handled by regexp).
% ------------------------- 0) SVG → disjoint cells -------------------------
LINE_WIDTH = 1.25;   % outline thickness for preview
NSEG       = 32;     % samples per Bezier segment (for Beziers)
AREA_EPS   = 1e-10;  % tiny area threshold
SIMPLIFY   = true;   % polyshape simplify

% --- read SVG
try
    xdoc = xmlread(svgFile);
catch ME
    error('Could not read SVG (%s): %s', svgFile, ME.message);
end

% --- collect closed loops from <polygon> and <path>
CellsXY = {};

% polygons
polys = xdoc.getElementsByTagName('polygon');
for i = 0:(polys.getLength-1)
    pts = char(polys.item(i).getAttribute('points'));
    if isempty(strtrim(pts)), continue; end
    nums = regexp(pts, '(-?\d*\.?\d+(?:[eE][-+]?\d+)?)', 'match');
    v    = str2double(nums);
    if numel(v) < 4, continue; end
    if mod(numel(v),2)~=0, v = v(1:end-1); end
    XY   = reshape(v,2,[])';
    if any(XY(1,:)~=XY(end,:)), XY = [XY; XY(1,:)]; end
    XY   = XY; %* PX2UM;                    % <<< scale to microns
    if size(XY,1) >= 4, CellsXY{end+1} = XY; end
end

% paths
paths = xdoc.getElementsByTagName('path');
for i = 0:(paths.getLength-1)
    d = char(paths.item(i).getAttribute('d'));
    if isempty(strtrim(d)), continue; end
    tokens = regexp(d,'([MmLlHhVvCcQqZz])|(-?\d*\.?\d+(?:[eE][-+]?\d+)?)','tokens');
    cmds = {}; vals = {};
    for t = 1:numel(tokens)
        pr = tokens{t};
        if ~isempty(pr{1}), cmds{end+1} = pr{1}; vals{end+1} = [];
        else, if isempty(vals), continue; end; vals{end} = [vals{end}, str2double(pr{2})];
        end
    end
    curr = [0 0]; startSub = [0 0]; pen = []; closed = {};
    for ci = 1:numel(cmds)
        c = cmds{ci}; v = vals{ci};
        switch c
            case {'M','m'}
                if mod(numel(v),2)~=0, v = v(1:end-1); end
                P = reshape(v,2,[])';
                for k = 1:size(P,1)
                    p = P(k,:); if c=='m', p = curr + p; end
                    if ~isempty(pen), pen = []; end
                    pen = p; startSub = p; curr = p;
                end
            case {'L','l'}
                if mod(numel(v),2)~=0, v = v(1:end-1); end
                P = reshape(v,2,[])';
                for k = 1:size(P,1), p = P(k,:); if c=='l', p = curr + p; end; pen=[pen; p]; curr=p; end
            case {'H','h'}
                for k=1:numel(v), x=v(k); if c=='h', x=curr(1)+x; end; p=[x,curr(2)]; pen=[pen; p]; curr=p; end
            case {'V','v'}
                for k=1:numel(v), y=v(k); if c=='v', y=curr(2)+y; end; p=[curr(1),y]; pen=[pen; p]; curr=p; end
            case {'C','c'}
                if mod(numel(v),6)~=0, v=v(1:6*floor(numel(v)/6)); end
                P0=curr;
                for k=1:6:numel(v)
                    cp1=[v(k) v(k+1)]; cp2=[v(k+2) v(k+3)]; P3=[v(k+4) v(k+5)];
                    if c=='c', cp1=P0+cp1; cp2=P0+cp2; P3=P0+P3; end
                    t  = linspace(0,1,max(2,round(NSEG)))';
                    w0=(1-t).^3; w1=3*(1-t).^2.*t; w2=3*(1-t).*t.^2; w3=t.^3;
                    B=[w0*P0(1)+w1*cp1(1)+w2*cp2(1)+w3*P3(1), ...
                       w0*P0(2)+w1*cp1(2)+w2*cp2(2)+w3*P3(2)];
                    if isempty(pen), pen=B; else, pen=[pen; B(2:end,:)]; end
                    P0=P3;
                end
                curr=P0;
            case {'Q','q'}
                if mod(numel(v),4)~=0, v=v(1:4*floor(numel(v)/4)); end
                P0=curr;
                for k=1:4:numel(v)
                    cp=[v(k) v(k+1)]; P2=[v(k+2) v(k+3)];
                    if c=='q', cp=P0+cp; P2=P0+P2; end
                    t=linspace(0,1,max(2,round(NSEG)))';
                    w0=(1-t).^2; w1=2*(1-t).*t; w2=t.^2;
                    B=[w0*P0(1)+w1*cp(1)+w2*P2(1), w0*P0(2)+w1*cp(2)+w2*P2(2)];
                    if isempty(pen), pen=B; else, pen=[pen; B(2:end,:)]; end
                    P0=P2;
                end
                curr=P0;
            case {'Z','z'}
                if ~isempty(pen)
                    if any(pen(1,:)~=pen(end,:)), pen = [pen; pen(1,:)]; end
                    if size(pen,1) >= 4
                        closed_subpaths{end+1} = pen;% * PX2UM;   % <<< scale to microns
                    end
                    pen = []; curr = startSub;
                end
        end
    end
    for k=1:numel(closed), CellsXY{end+1}=closed{k}; end %#ok<SAGROW>
end
fprintf('Found %d closed loops (candidate cells).\n', numel(CellsXY));

%% ------------------ 3) CALIBRATION: SVG UNITS → MICRONS -------------------
% Purpose:
%   - Shift coordinates by viewBox origin so (0,0) = SVG canvas top-left.
%   - Apply per-axis scaling PX2UM_X, PX2UM_Y to get µm.
% Why:
%   - Downstream areas/thickness should be physically meaningful.
% Pitfalls:
%   - If no viewBox: minX=minY=0 assumed.
%   - Anisotropy: use different PX2UM_X vs PX2UM_Y if needed.
%
% === Calibration: user-units → microns ===
% If the SVG defines a viewBox, shift by its origin (minX,minY) so (0,0) = canvas top-left
svgRoot = xdoc.getDocumentElement();
vbAttr  = char(svgRoot.getAttribute('viewBox'));   % "minX minY width height"
if ~isempty(strtrim(vbAttr))
    vb   = sscanf(vbAttr, '%f %f %f %f');
    minX = vb(1);  minY = vb(2);
else
    % No viewBox → assume already in canvas coordinates
    minX = 0;      minY = 0;
end

% Apply shift + scale to every loop
for i = 1:numel(CellsXY)
    CellsXY{i}(:,1) = (CellsXY{i}(:,1) - minX) * PX2UM_X;
    CellsXY{i}(:,2) = (CellsXY{i}(:,2) - minY) * PX2UM_Y;
end

fprintf('Applied constant scale SX=%.6f, SY=%.6f µm/unit (viewBox origin: %.3f, %.3f)\n', ...
        PX2UM_X, PX2UM_Y, minX, minY);

%% ------------- 4) DISJOINT POLYSHAPES (OVERLAP RESOLUTION) ----------------
% Purpose:
%   - Clean each loop into polyshape regions; drop tiny/invalid areas.
%   - Sort by area (largest→smallest) and greedily subtract the union so cells don't overlap.
% Why:
%   - A consistent, non-overlapping cell set ('Pcells') is needed for meshing & stats.
% Knobs:
%   - SIMPLIFY (polyshape's internal simplifier), AREA_EPS (tiny area threshold).
% Diagnostics:
%   - Prints cell count kept after resolution.
%
% --- make disjoint polyshapes
Ptmp  = polyshape.empty(0,1); areas = [];
for i=1:numel(CellsXY)
    XY=CellsXY{i};
    if size(XY,1)>=2 && all(XY(1,:)==XY(end,:)), XY=XY(1:end-1,:); end
    if size(XY,1)<3, continue; end
    P = polyshape(XY(:,1), XY(:,2), 'Simplify', SIMPLIFY);
    if P.NumRegions==0 || area(P)<=AREA_EPS, continue; end
    regs=regions(P);
    for r=1:numel(regs)
        ar=area(regs(r));
        if ar>AREA_EPS, Ptmp(end+1,1)=regs(r); areas(end+1,1)=ar; end %#ok<SAGROW>
    end
end
if isempty(Ptmp), error('No valid cell polygons created.'); end
[~,ord]=sort(areas,'descend'); Ptmp=Ptmp(ord);

Pcells = polyshape.empty(0,1);
U = polyshape();
for i=1:numel(Ptmp)
    Pi = subtract(Ptmp(i), U);
    if area(Pi)>AREA_EPS
        regs=regions(Pi);
        for r=1:numel(regs)
            if area(regs(r))>AREA_EPS
                Pcells(end+1,1)=regs(r); U=union(U,regs(r)); %#ok<SAGROW>
            end
        end
    end
end
Nc = numel(Pcells);
fprintf('Kept %d disjoint cells after overlap resolution.\n', Nc);


%% ------ 5) BANDS AROUND WALLS + GLOBAL MESH (SHARED NODES) ----------------
% Purpose:
%   - Build per-cell 'wall bands' (ring = cell \ buffer(cell, -w/2)).
%   - Snap all band boundary points to a global grid (TOL) so adjacent cells share nodes.
%   - Constrained Delaunay triangulation over all bands at once.
%   - Classify boundary edges as 'inner' (faces lumen) vs 'exterior' (faces outside tissue).
% Why:
%   - Shared nodes enforce continuity along shared walls (no cracks).
%   - Edge classes let you apply different loads/props on outer rim vs inner walls.
% Knobs:
%   - w (target in-plane band width) auto-set from median cell size; auto-thins on tiny cells.
%   - hmax (boundary densify spacing); smaller → finer mesh.
%   - TOL (node snapping tolerance) ~ 1e-5 * typLen is a good start.
% Pitfalls:
%   - If band creation fails for small cells, SHRINK_TRIES / SHRINK_FACTOR will thin w.
%   - Classification uses a small probe along outward normal; epsn should be ≤ ~w/4.
% Quick preview:
%   - After classification, plot exterior edges in red to visually validate tissue perimeter.

% ------------------------- 1) Bands per cell & meshing -------------------------
% --- 1a) Bands just to get w_per_cell (reuse your loop) ---
% Result: w_per_cell (Nc x 1). You can set a fallback if you prefer:
% w_per_cell = 0.03*sqrt(median(arrayfun(@area,Pcells))) * ones(Nc,1);
cellAreas = arrayfun(@area, Pcells);
typLen    = sqrt(median(cellAreas));
w         = 0.03*typLen;          % in-plane band thickness (global target)
hmax      = 0.5*w;                % boundary point spacing for meshing
TOL       = max(1e-5*typLen,1e-9);
SHRINK_TRIES  = 6;
SHRINK_FACTOR = 0.7;

% 1a) build interior bands Bi = cell \ buffer(cell,-w/2), with auto-thin for tiny cells
B = polyshape.empty(0,1);
w_per_cell = w*ones(Nc,1);
for i=1:Nc
    wi = w/2; Bi = polyshape(); ok=false;
    for attempt=1:SHRINK_TRIES
        Pin = polybuffer(Pcells(i), -wi);
        if Pin.NumRegions>0 && area(Pcells(i))-area(Pin) > 1e-12
            Bi = subtract(Pcells(i), Pin); ok=true; break
        else
            wi = wi * SHRINK_FACTOR;
        end
    end
    if ~ok, error('Band creation failed for cell %d. Reduce w.', i); end
    w_per_cell(i) = 2*wi; B(i,1)=Bi;
end
fprintf('Bands built. Median thickness = %.4g (min %.4g, max %.4g)\n', ...
        median(w_per_cell), min(w_per_cell), max(w_per_cell));

% --- Build Lumen polygons from Pcells + w_per_cell ---
if ~exist('Lumen','var') || isempty(Lumen)
    Lumen = polyshape.empty(0,1);
    for ci = 1:Nc
        Lumen(ci,1) = polybuffer(Pcells(ci), -0.5*w_per_cell(ci));
    end
end

% 1b) collect all band boundaries → global XY, CE (with snapping to share nodes)
key2idx = containers.Map('KeyType','char','ValueType','double');
addPt = @(p) add_point_snap(p, TOL, key2idx);
XY = []; CE = [];
for ci=1:Nc
    [x,y] = boundary(B(ci)); if isempty(x), continue; end
    isn = isnan(x)|isnan(y); cuts=[0; find(isn); numel(x)+1];
    for k=1:numel(cuts)-1
        a=cuts(k)+1; b=cuts(k+1)-1; if b<=a, continue; end
        loop=[x(a:b) y(a:b)]; if ~isequal(loop(1,:),loop(end,:)), loop=[loop; loop(1,:)]; end
        loop = densifyPolyline(loop, hmax);
        idx = zeros(size(loop,1),1);
        for m=1:size(loop,1), idx(m)=addPt(loop(m,:)); end
        CE = [CE; [idx(1:end-1) idx(2:end)]]; %#ok<SAGROW>
    end
end
XY = materialize_map_points(key2idx);
CE = unique(sort(CE,2),'rows');

% 1c) single constrained triangulation over all bands
dt = delaunayTriangulation(XY, CE);
T  = dt.ConnectivityList;
P  = dt.Points;

% keep triangles whose centroid is inside any band; tag with band_id
ctr = (P(T(:,1),:)+P(T(:,2),:)+P(T(:,3),:))/3;
keep = false(size(T,1),1); band_id = zeros(size(T,1),1);
for ci=1:Nc
    IN = isinterior(B(ci), ctr(:,1), ctr(:,2));
    band_id(~keep & IN) = ci;
    keep = keep | IN;
end
T = T(keep,:); band_id = band_id(keep);
fprintf('Band mesh: %d nodes, %d triangles.\n', size(P,1), size(T,1));

% 1d) boundary edges of the band mesh + outward normals + classify inner/exterior
[Eall,E2T] = tri_edges_attach(T);
isBdry = (E2T(:,2)==0);
Ebdry  = Eall(isBdry,:);        % boundary edges [i j]
Tab    = E2T(isBdry,1);         % adjacent triangle (inside the band)
nBdry  = size(Ebdry,1);

% outward normal (away from triangle interior)
n_out = zeros(nBdry,2);
for r=1:nBdry
    i=Ebdry(r,1); j=Ebdry(r,2); t=Tab(r);
    pij = P(j,:)-P(i,:); L=hypot(pij(1),pij(2)); if L==0, continue; end
    tvec = pij / L; nL = [-tvec(2), tvec(1)];
    tri=T(t,:); k3=tri(~ismember(tri,[i j]));
    s = tvec(1)*(P(k3,2)-P(i,2)) - tvec(2)*(P(k3,1)-P(i,1));
    n_out(r,:) = (s>0) * (-nL) + (s<=0) * (nL);
end

% classify boundary edges: inner (faces the same cell's lumen) vs exterior
epsn = 0.25*median(w_per_cell);
Uall = polyshape(); for ci=1:Nc, Uall = union(Uall, Pcells(ci)); end
edgeCell  = zeros(nBdry,1);
isInner   = false(nBdry,1);
isExterior= false(nBdry,1);

for r=1:nBdry
    t=Tab(r); c=band_id(t); edgeCell(r)=c;
    i=Ebdry(r,1); j=Ebdry(r,2); mid=0.5*(P(i,:)+P(j,:));
    probe = mid + epsn*n_out(r,:);
    if isinterior(Pcells(c), probe(1), probe(2))
        isInner(r)=true;       % lumen-side edge of cell c
    elseif ~isinterior(Uall, probe(1), probe(2))
        isExterior(r)=true;    % outer perimeter of tissue
    end
end

%% -------- 6) OPTIONAL: OUT-OF-PLANE THICKNESS MAPS (t_elem, t_edge) -------
% Purpose:
%   - Define per-element (t_elem) and per-boundary-edge (t_edge) thickness multipliers.
%   - Example: bump outer perimeter triangles/edges.
% Why:
%   - Useful when keeping geometry fixed but adjusting stiffness/loads by weights.
%
% ------------------------- 2) Optional Variable (out of plane) thickness maps -------------------------
% base thickness for all walls; outer epidermis walls = 8× base
t_base  = 1.0;
t_outer = 1.0 * t_base;

nTri  = size(T,1);

% per-triangle thickness for stiffness: triangles adjacent to exterior edges → 8×
t_elem = t_base * ones(nTri,1);
t_elem(unique(Tab(isExterior))) = t_outer;  % any tri touching exterior edge is thick

% per-boundary-edge thickness for line loads (if you ever use p_ext)
t_edge = t_base * ones(nBdry,1);
t_edge(isExterior) = t_outer;


%% ------------------------- VIS: who gets 2× thickness? -------------------------
% Triangles using 8× thickness = any tri adjacent to an exterior edge
tri_thick = unique(Tab(isExterior));
tri_thick = tri_thick(tri_thick > 0);

% Cells (bands) affected = band_id of those triangles
cells_epid = unique(band_id(tri_thick));

fprintf('Epidermis (2×) is applied to %d cells.\n', numel(cells_epid));
% Uncomment to print the actual indices:
% disp('Cell indices with 2× thickness:'); disp(cells_epid.');

% % (A) Mesh-level view: orange triangles are 8×; red edges are exterior
% figure('Color','w'); hold on; axis equal off
% triplot(T, P(:,1), P(:,2), 'Color', [0.85 0.85 0.85]);                  % all bands (light)
% if ~isempty(tri_thick)
%     patch('Faces',T(tri_thick,:),'Vertices',P, ...
%           'FaceColor',[1.0 0.80 0.50],'EdgeColor','none','FaceAlpha',0.9);
% end
% for r = find(isExterior).'
%     e = Ebdry(r,:);
%     plot(P(e,1), P(e,2), 'r-', 'LineWidth', 1.5);                       % outer edges
% end
% title('2× thickness: orange = thick triangles, red = exterior edges')

% % (B) Cell-level view: blue outlines mark cells that have any 8× triangles
% figure('Color','w'); hold on; axis equal off
% % draw all cells (thin gray)
% for ci = 1:numel(Pcells)
%     plot(Pcells(ci),'FaceColor','none','EdgeColor',[0.7 0.7 0.7],'LineWidth',0.5);
% end
% % highlight epidermis cells (bold blue) and label their indices
% for ci = cells_epid.'
%     plot(Pcells(ci),'FaceColor','none','EdgeColor',[0 0.40 0.80],'LineWidth',1.6);
%     [cx,cy] = centroid(Pcells(ci));
%     text(cx,cy, sprintf('%d',ci), 'Color',[0 0.30 0.60], ...
%          'FontSize',8,'HorizontalAlignment','center');
% end
% title(sprintf('Cells receiving 2× thickness (count = %d)', numel(cells_epid)))

% (Optional) heat-map-ish summary: fraction of each cell''s band that is 8×
tris_per_cell       = accumarray(band_id, 1, [Nc,1], @sum, 0);
tris_thick_per_cell = accumarray(band_id(tri_thick), 1, [Nc,1], @sum, 0);
frac_thick = tris_thick_per_cell ./ max(1, tris_per_cell);  % in [0,1]

% figure('Color','w'); hold on; axis equal off
% for ci = 1:Nc
%     f = frac_thick(ci);
%     if f==0
%         fc = [1 1 1];     % white = none
%     else
%         fc = [1, 1-0.8*f, 1-0.8*f];  % pale→deep salmon as coverage increases
%     end
%     plot(Pcells(ci),'FaceColor',fc,'EdgeColor',[0.2 0.2 0.2],'LineWidth',0.5,'FaceAlpha',0.7);
% end
% title('Per-cell coverage of 2× triangles (white = none, deeper = more)')
% colorbar('off'); % purely categorical shading here

%% ---- 7) RADIAL (IN-PLANE) THICK EPIDERMIS + RE-MESH (GEOMETRY-ENCODED) ---
% Purpose:
%   - Construct a global outer belt (width w_out) inside the tissue and union it into
%     only the peripheral portions of each cell band.
%   - Rebuild the global constrained triangulation on the widened bands (B2).
%   - Recompute boundary classification; rebuild exact lumens as Pcells \ B2.
% Why:
%   - Encodes thickness *geometrically* (constant FE thickness t0 thereafter).
% Knobs:
%   - xThickness: multiplies base band width (w_in) to set the outer rim width (w_out).
%   - Auto-thinning tries if belt can't be built in very narrow regions.
% Pitfalls:
%   - Very thin epidermal arcs may require lowering xThickness or increasing SHRINK_TRIES.
% Convention:
%   - Once geometry encodes thickness, set FE thickness t0 constant and DO NOT use t_elem/t_edge multipliers.
%
% ------------------------- 2) Radial (in-plane) thick epidermis + re-mesh -------------------------
% Goal: make the outer wall 8× thicker *in-plane* (radially inward from tissue perimeter).
% Keeps a CONSTANT FE thickness later (no t_elem/t_edge multipliers).

xThickness = 2;

% Base inner-wall thickness (use what you actually built for the bands)
if exist('w_per_cell','var') && ~isempty(w_per_cell)
    w_in = median(w_per_cell);         % your current band thickness ~ in-plane width
else
    w_in = 0.03*typLen;                % fallback if bands were not built that way
end
w_out = xThickness * w_in;                      % desired radial thickness on the outer rim

% (a) Outer belt of width w_out inside the tissue
Uall = polyshape(); for ci=1:Nc, Uall = union(Uall, Pcells(ci)); end
SHRINK_TRIES  = 6; SHRINK_FACTOR = 0.7;
belt_ok = false; wtry = w_out;
for attempt = 1:SHRINK_TRIES
    Uin = polybuffer(Uall, -wtry);                 % shrink tissue inward by wtry
    if Uin.NumRegions > 0 && area(Uall) - area(Uin) > 1e-12
        G_outer = subtract(Uall, Uin);             % the interior belt (width ≈ w_out)
        belt_ok = true; break
    else
        wtry = wtry * SHRINK_FACTOR;               % auto-thin if geometry is too narrow
    end
end
if ~belt_ok
    warning('Could not build the outer belt at w_out=%.3g; skipping radial boost.', w_out);
    G_outer = polyshape();   % noop
end

% (b) Update each cell's band: add the part of the global belt that lies in that cell
%     This widens only the epidermal arcs; interior cell–cell walls remain at base width.
B2 = B;              % start from your existing per-cell bands (base ~ w_in)
epi_idx = [];        % record which cells actually got widened
for ci = 1:Nc
    extra = intersect(G_outer, Pcells(ci));        % radial-thick strip inside this cell
    if extra.NumRegions>0 && area(extra) > 0
        B2(ci) = union(B2(ci), extra);
        epi_idx(end+1) = ci; %#ok<AGROW>
    end
end
fprintf('Radial thickening applied. Base=%.3g, outer=%.3g. Epidermis cells auto-detected: %d\n', ...
        w_in, w_out, numel(epi_idx));

% (c) Re-collect ALL band boundaries from B2 → global node list & constraints, then re-mesh
key2idx = containers.Map('KeyType','char','ValueType','double');
addPt = @(p) add_point_snap(p, TOL, key2idx);
XY = []; CE = [];

for ci = 1:Nc
    [x,y] = boundary(B2(ci));
    if isempty(x), continue; end
    isn = isnan(x) | isnan(y);
    cuts = [0; find(isn); numel(x)+1];
    for k = 1:numel(cuts)-1
        a = cuts(k)+1; b = cuts(k+1)-1;
        if b<=a, continue; end
        loop = [x(a:b) y(a:b)];
        if ~isequal(loop(1,:), loop(end,:)), loop = [loop; loop(1,:)]; end
        loop = densifyPolyline(loop, hmax);
        idx = zeros(size(loop,1),1);
        for m = 1:size(loop,1), idx(m) = addPt(loop(m,:)); end
        CE  = [CE; [idx(1:end-1) idx(2:end)]]; %#ok<AGROW>
    end
end

% Materialize nodes and unique constraints
P = materialize_map_points(key2idx);
CE = unique(sort(CE,2),'rows');

% Constrained Delaunay triangulation
dt = delaunayTriangulation(P, CE);
T  = dt.ConnectivityList;
P  = dt.Points;

% Keep only triangles inside any band; tag owner band_id
ctr = (P(T(:,1),:)+P(T(:,2),:)+P(T(:,3),:))/3;
keep = false(size(T,1),1); band_id = zeros(size(T,1),1);
for ci = 1:Nc
    IN = isinterior(B2(ci), ctr(:,1), ctr(:,2));
    band_id(~keep & IN) = ci;
    keep = keep | IN;
end
T = T(keep,:); band_id = band_id(keep);
fprintf('Re-meshed (radially thickened) bands: %d nodes, %d triangles.\n', size(P,1), size(T,1));

% (d) Boundary edges and classification (inner vs exterior) on the new mesh
[Eall, E2T] = tri_edges_attach(T);
isBdry = (E2T(:,2)==0);
Ebdry  = Eall(isBdry,:);         % boundary edges [i j]
Tab    = E2T(isBdry,1);          % adjacent tri index
nBdry  = size(Ebdry,1);

% Outward normals (from band interior to outside)
n_out = zeros(nBdry,2);
for r = 1:nBdry
    i = Ebdry(r,1); j = Ebdry(r,2); t = Tab(r);
    pij = P(j,:)-P(i,:); L = hypot(pij(1),pij(2)); if L==0, continue; end
    tvec = pij / L; nL = [-tvec(2), tvec(1)];
    tri = T(t,:); k3 = tri(~ismember(tri,[i j]));
    s = tvec(1)*(P(k3,2)-P(i,2)) - tvec(2)*(P(k3,1)-P(i,1));
    n_out(r,:) = (s>0) * (-nL) + (s<=0) * (nL);
end

% Classify boundary edges
epsn = 0.25 * w_in;               % small probe inward
Uall = polyshape(); for ci=1:Nc, Uall = union(Uall, Pcells(ci)); end
edgeCell  = zeros(nBdry,1);
isInner   = false(nBdry,1);
isExterior= false(nBdry,1);
for r = 1:nBdry
    t = Tab(r); c = band_id(t); edgeCell(r) = c;
    i = Ebdry(r,1); j = Ebdry(r,2); mid = 0.5*(P(i,:)+P(j,:));
    probe = mid + epsn*n_out(r,:);
    if isinterior(Pcells(c), probe(1), probe(2))
        isInner(r) = true;                 % lumen-side
    elseif ~isinterior(Uall, probe(1), probe(2))
        isExterior(r) = true;              % outer perimeter
    end
end

% (e) Rebuild Lumen from geometry (exact): lumen = cell \ widened band
Lumen = polyshape.empty(0,1);
for ci = 1:Nc
    Lumen(ci,1) = subtract(Pcells(ci), B2(ci));
end

% (f) FE thickness now CONSTANT (geometry encodes radial thickening)
t0 = 1.0;                % use t0 in stiffness & tractions later
% (Downstream changes:)
%   - In element stiffness:  ke = t0 * Atri * (Bm.' * D * Bm);
%   - In inner-edge traction: fnode = (t0 * Lij / 2) * q;
%   - Do NOT use t_elem / t_edge multipliers anymore.

%% ------------ 8) VIS SANITY: ADDED RADIAL STRIPS & TISSUE BOUNDARY --------
% Purpose:
%   - Show only the *added* outer strips (B2 \ B) on top of cell outlines.
%   - Label affected cells; draw tissue boundary for quick eyeballing of SVG quality.
% Why:
%   - Fast "does my SVG need cleanup?" visual: gaps or stray edges stand out immediately.
%  --- VIS: Added radial strips only (epidermal widening sanity check) ---
% Shows B2 \ B (i.e., the extra belt that made the outer rim thicker).
% Falls back to intersect(G_outer, Pcells(ci)) if B is not available.

% guards / fallbacks
if ~exist('Uall','var') || isempty(Uall)
    Uall = polyshape(); for ci=1:Nc, Uall = union(Uall, Pcells(ci)); end
end
if ~exist('w_in','var'),  w_in  = 0.03*sqrt(median(arrayfun(@area,Pcells))); end
if ~exist('w_out','var'), w_out = xThickness*w_in; end

figure('Color','w'); hold on; axis equal off
ax = gca; set(ax,'YDir','reverse'); 
title(sprintf('Added radial strips only  (base=%.3g, outer=%.3g)', w_in, w_out));

% light outlines of all cells
for ci = 1:Nc
    plot(Pcells(ci), 'FaceColor','none', 'EdgeColor',[0.75 0.75 0.75], 'LineWidth',0.5);
end

colAdd = [1.00 0.70 0.20];  % orange for added strip
nAdd = 0;

for ci = 1:Nc
    % "added" = widened part of the band inside this cell
    if exist('B','var') && ci <= numel(B) && B(ci).NumRegions>0
        added = subtract(B2(ci), B(ci));         % preferred: exact delta
    else
        if exist('G_outer','var') && ~isempty(G_outer)
            added = intersect(G_outer, Pcells(ci)); % fallback: belt ∩ cell
        else
            added = polyshape(); % nothing to show
        end
    end

    if added.NumRegions>0 && area(added) > 0
        plot(added, 'FaceColor', colAdd, 'EdgeColor','none', 'FaceAlpha', 0.9);
        nAdd = nAdd + 1;

        % label the cell index at its centroid (helps double-check membership)
        [cx,cy] = centroid(Pcells(ci));
        text(cx,cy,num2str(ci), 'Color',[0.25 0.25 0.25], 'FontSize',8, ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontWeight','bold','Clipping','on');
    end
end

% tissue boundary
plot(Uall, 'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1.0);
text(double(mean(xlim)), double(max(ylim)), sprintf('N_{epi with added strip} = %d', nAdd), ...
     'HorizontalAlignment','center','VerticalAlignment','top','Color',[0.2 0.2 0.2]);

% optional legend anchors
plot(NaN,NaN,'-','Color',[0.75 0.75 0.75],'LineWidth',0.5);  % cell outlines
patch(NaN,NaN,[1 0.7 0.2], 'EdgeColor','none');              % added strip
plot(NaN,NaN,'-k','LineWidth',1.0);                          % tissue boundary
legend({'cells','added strip','tissue boundary'}, 'Location','bestoutside');

%% =============== 9) INTERACTIVE: CLICK LUMENS TO SET CELL TYPES ===========
% Purpose:
%   - Click lumens to assign epidermis → cortex → middle cortex → endodermis; remainder = stele.
% Why:
%   - Downstream area stats and FE can depend on cell type.
% UI:
%   - Left-click add, right-click remove, Enter/Esc to finish each group.
% Fallback:
%   - If a lumen is empty for a cell, falls back to whole cell for visualization only.
% ========================= Interactive: click lumens to assign cell types =========================
% Returns:
%   epidermis, cortex, middle cortex, endodermis, stele   (index vectors, 1..Nc)

if ~exist('Nc','var'), Nc = numel(Pcells); end
if ~exist('Lumen','var') || isempty(Lumen)
    % build lumens just inside the wall bands
    Lumen = polyshape.empty(0,1);
    for ci = 1:Nc
        Lumen(ci,1) = polybuffer(Pcells(ci), -0.5*w_per_cell(ci));
    end
end

[epidermis, cortex, middle_cortex, endodermis, stele] = click_select_cell_types(Lumen, Pcells);

% print nicely
fprintf('\nSelected groups:\n');
fprintf('  epidermis  = %s\n', mat2str(epidermis));
fprintf('  cortex     = %s\n', mat2str(cortex));
fprintf('  middle_cortex     = %s\n', mat2str(middle_cortex));
fprintf('  endodermis = %s\n', mat2str(endodermis));
fprintf('  stele      = %s\n\n', mat2str(stele));

% ------------------------------------------------------------------------
function [epi, cor, midc, endo, stele] = click_select_cell_types(Lumen, Pcells)
    Nc = numel(Pcells);
    % colors
    colBG = [0.96 0.96 0.96]; colTxt = [0.2 0.2 0.2];
    colE  = [0.10 0.45 0.85]; % epidermis (blue)
    colC  = [0.10 0.65 0.10]; % cortex    (green)
    colM  = [0.55 0.10 0.70]; % middle cortex (purple)
    colN  = [0.70 0.40 0.05]; % endodermis(brown)
    colS  = [0.70 0.70 0.70]; % unassigned/stele (gray)

    % base figure with all lumens and index labels
    fig = figure('Color','w','Name','Click-to-assign cell types','NumberTitle','off');
    ax  = axes('Parent',fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
    title(ax, 'Click lumens to assign: Epidermis → Cortex → Middle Cortex → Endodermis.  (Enter to finish each)', ...
          'FontWeight','bold');

    % fills and outlines
    hFill = gobjects(Nc,1);
    for ci = 1:Nc
        if Lumen(ci).NumRegions>0
            hFill(ci) = plot(Lumen(ci), 'FaceColor', colBG, 'EdgeColor','none', 'FaceAlpha', 0.8, 'Parent', ax);
        else
            % fallback to whole cell if lumen missing
            hFill(ci) = plot(Pcells(ci), 'FaceColor', colBG, 'EdgeColor','none', 'FaceAlpha', 0.8, 'Parent', ax);
        end
        plot(Pcells(ci), 'FaceColor','none', 'EdgeColor',[0.3 0.3 0.3], 'LineWidth', 0.6, 'Parent', ax);
    end

    % numeric cell indices at lumen centroids
    for ci = 1:Nc
        if Lumen(ci).NumRegions>0
            [cx,cy] = centroid(Lumen(ci));
        else
            [cx,cy] = centroid(Pcells(ci));
        end
        text(cx, cy, num2str(ci), 'Color', colTxt, 'FontSize', 8, ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontWeight','bold', 'Clipping','on', 'Parent', ax);
    end

    % assign per type (order matters)
    taken = false(Nc,1);
    epi  = select_one_type(ax, Lumen, 'EPIDERMIS',     colE, taken, hFill);
    taken(epi) = true;

    cor  = select_one_type(ax, Lumen, 'CORTEX',        colC, taken, hFill);
    taken(cor) = true;

    midc = select_one_type(ax, Lumen, 'MIDDLE CORTEX', colM, taken, hFill);
    taken(midc) = true;

    endo = select_one_type(ax, Lumen, 'ENDODERMIS',    colN, taken, hFill);
    taken(endo) = true;

    % the rest → stele
    stele = find(~taken).';
    % color remaining as stele
    set(hFill(stele), {'FaceColor'}, num2cell(repmat(colS, numel(stele),1),2));

    % final legend-ish title
    title(ax, sprintf('Done. Epi:%d  Cort:%d  MidC:%d  Endo:%d  Stele:%d', ...
                      numel(epi), numel(cor), numel(midc), numel(endo), numel(stele)));
end

function idx = select_one_type(ax, Lumen, typeName, col, taken, hFill)
    % Interactive selection loop for one type.
    % Left-click = add; Right-click = remove; Enter/Esc = finish.
    Nc = numel(Lumen);
    idx = []; sel = false(Nc,1);

    instr = sprintf('Select %s:  Left-click=add, Right-click=remove, Enter/Esc=finish', typeName);
    title(ax, instr, 'FontWeight','bold', 'Color', col);

    while ishghandle(ax)
        try
            [x,y,button] = ginput(1);    % mouse or key
        catch
            break;                       % figure closed
        end
        if isempty(button) || any(button==[13 27])    % Enter(13) or Esc(27)
            break;
        end

        % find which lumen contains the point
        hit = 0;
        for ci = 1:Nc
            if Lumen(ci).NumRegions>0 && isinterior(Lumen(ci), x, y)
                hit = ci; break;
            end
        end
        if hit==0
            % clicked empty space
            continue;
        end
        if taken(hit) && ~sel(hit)
            % already assigned to a previous type
            % quick flash to indicate blocked
            old = hFill(hit).FaceColor; set(hFill(hit),'FaceColor',[1 0.6 0.6]); drawnow; pause(0.05);
            set(hFill(hit),'FaceColor',old); drawnow;
            continue;
        end

        switch button
            case 1   % left: add/select
                sel(hit) = true;
                set(hFill(hit), 'FaceColor', col, 'FaceAlpha', 0.75);
            case 3   % right: remove/unselect
                sel(hit) = false;
                set(hFill(hit), 'FaceColor', [0.96 0.96 0.96], 'FaceAlpha', 0.8);
        end

        idx = find(sel).';
        title(ax, sprintf('%s selected: %d   (Enter/Esc to finish)', typeName, numel(idx)), 'Color', col);
        drawnow;
    end
end


%% ------ 10) PER-TYPE LUMEN AREAS & DEVIATION FROM TYPE MEAN (A0) ----------
% Purpose:
%   - Compute lumen area per cell (A0), type-wise means, and each cell's % deviation.
% Why:
%   - Quantify size scaling across layers, or identify outliers within a type.
% Output:
%   - Console printouts of counts and means; 'dev_pct' per cell for plotting/exports.
% ------------------------- 4d) Cell-type selection + area deviation from type mean -------------------------
% Build a type map: 1=Epidermis, 2=Cortex, 3=Endodermis, 4=Stele
cell_type = zeros(Nc,1);
cell_type(epidermis)     = 1;
cell_type(cortex)        = 2;
cell_type(middle_cortex) = 3;
cell_type(endodermis)    = 4;
cell_type(stele)         = 5;
type_names = {'Epidermis','Cortex','Middle Cortex','Endodermis','Stele'};


% Ensure we have lumen polygons (undeformed) to measure area
if ~exist('Lumen','var') || isempty(Lumen)
    Lumen = polyshape.empty(0,1);
    for ci = 1:Nc
        Lumen(ci,1) = polybuffer(Pcells(ci), -0.5*w_per_cell(ci));
    end
end

% Per-cell lumen area (undeformed)
A0 = NaN(Nc,1);
for ci = 1:Nc
    if Lumen(ci).NumRegions > 0
        A0(ci) = area(Lumen(ci));
    end
end

% Per-type mean area (ignore NaNs / empty lumens)
T = 5;  % number of types
type_mean  = NaN(T,1);
type_count = zeros(T,1);
for t = 1:T
    idx = find(cell_type==t & ~isnan(A0) & A0>0);
    type_count(t) = numel(idx);
    if ~isempty(idx), type_mean(t) = mean(A0(idx)); end
end

% Per-cell signed percent deviation from its type mean: 100*(Ai - Atype)/Atype
dev_pct = NaN(Nc,1);
for ci = 1:Nc
    t = cell_type(ci);
    if t>0 && ~isnan(type_mean(t)) && type_mean(t) > 0 && ~isnan(A0(ci))
        dev_pct(ci) = 100 * (A0(ci) - type_mean(t)) / type_mean(t);
    end
end

% Console summary
fprintf('\nAverage lumen area per type:\n');
for t = 1:4
    fprintf('  %-11s  N=%3d   mean(A0)=%.4g\n', type_names{t}, type_count(t), type_mean(t));
end
fprintf('Global lumen area (all cells): mean=%.4g, N=%d\n\n', mean(A0(~isnan(A0))), sum(~isnan(A0)));

%%
%% ---- Lumen area bar plot (all types) with error bars ----
% Requires indices: epidermis, cortex, middle_cortex, endodermis, stele
% Also needs: Pcells, Nc, and (optionally) Lumen & w_per_cell
% Error bars: set err_mode = 'sem' or 'std'

err_mode = 'sem';   % 'sem' (default) or 'std'
ylab     = 'Lumen area (model units^2)';  % change if you know it's µm^2

if ~exist('Nc','var'), Nc = numel(Pcells); end

% Build lumens if not already available
if ~exist('Lumen','var') || isempty(Lumen)
    if ~exist('w_per_cell','var') || isempty(w_per_cell)
        % fallback thickness guess if bands weren't built
        cellAreas = arrayfun(@area, Pcells);
        typLen    = sqrt(median(cellAreas));
        w_per_cell = 0.03*typLen * ones(Nc,1);
    end
    Lumen = polyshape.empty(0,1);
    for ci = 1:Nc
        Lumen(ci,1) = polybuffer(Pcells(ci), -0.5*w_per_cell(ci));
    end
end

% Check groups exist (define empty if missing so code still runs)
if ~exist('epidermis','var'),     epidermis = []; end
if ~exist('cortex','var'),        cortex = []; end
if ~exist('middle_cortex','var'), middle_cortex = []; end
if ~exist('endodermis','var'),    endodermis = []; end
if ~exist('stele','var'),         stele = []; end

type_names = {'Epidermis','Cortex','Middle Cortex','Endodermis','Stele'};
groups     = {epidermis, cortex, middle_cortex, endodermis, stele};

% Colors (match your interactive picker palette)
colE = [0.10 0.45 0.85]; % epidermis (blue)
colC = [0.10 0.65 0.10]; % cortex    (green)
colM = [0.00 0.60 0.55]; % middle cortex (teal)
colN = [0.70 0.40 0.05]; % endodermis (brown)
colS = [0.70 0.70 0.70]; % stele (gray)
cols = [colE; colC; colM; colN; colS];

% Helper to collect areas for a group
getAreas = @(idx) collect_areas(idx, Lumen, Nc);

areas    = cellfun(getAreas, groups, 'UniformOutput', false);
n_vec    = cellfun(@(a) sum(isfinite(a) & a>0), areas);
mu_vec   = cellfun(@(a) mean(a, 'omitnan'), areas);
sd_vec   = cellfun(@(a) std(a,  'omitnan'), areas);
sem_vec  = sd_vec ./ sqrt(max(n_vec,1));
err_vec  = strcmpi(err_mode,'std') .* sd_vec + strcmpi(err_mode,'sem') .* sem_vec;

% Replace NaNs for empty groups to plot cleanly
mu_plot  = mu_vec;  mu_plot(~isfinite(mu_plot)) = 0;
err_plot = err_vec; err_plot(~isfinite(err_plot)) = 0;

% --- Bar plot with error bars ---
figure('Color','k'); hold on
b = bar(mu_plot, 'FaceColor', 'flat');
b.CData = cols;   % color each bar by type
er = errorbar(1:numel(mu_plot), mu_plot, err_plot, 'k', 'LineStyle','none', 'LineWidth', 1.2);

% Optional: overlay individual data points (jittered)
do_points = true;
if do_points
    for t = 1:numel(areas)
        ai = areas{t};
        ai = ai(isfinite(ai) & ai>0);
        if isempty(ai), continue; end
        xj = t + 0.06*randn(size(ai));
        scatter(xj, ai, 18, 'filled', 'MarkerFaceColor', cols(t,:), ...
            'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','none');
    end
end

% Annotate N above bars
for t = 1:numel(mu_plot)
    if n_vec(t) > 0
        ytxt = mu_plot(t) + (err_plot(t) + 1e-12)*1.08;
        text(t, ytxt, sprintf('N=%d', n_vec(t)), ...
             'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
             'FontSize', 8, 'Color',[0.2 0.2 0.2]);
    else
        % mark empty group subtly at y=0
        text(t, 0, 'N=0', 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontSize',8, 'Color',[0.5 0.5 0.5]);
    end
end

% Axes cosmetics
xlim([0.5, numel(mu_plot)+0.5]);
set(gca, 'XTick', 1:numel(mu_plot), 'XTickLabel', type_names, 'XTickLabelRotation', 20);
ylabel(ylab);
title(sprintf('Lumen area by type (bars = mean, error = %s)', upper(err_mode)));
grid on; box on

% Console summary
fprintf('\nLumen area by type (%s):\n', upper(err_mode));
for t = 1:numel(type_names)
    fprintf('  %-13s  N=%3d  mean=%.4g  std=%.4g  sem=%.4g\n', ...
        type_names{t}, n_vec(t), mu_vec(t), sd_vec(t), sem_vec(t));
end

%% ------------------------ 11) HELPER FUNCTIONS ----------------------------
% densifyPolyline  : subdivide boundary segments to cap max edge length ~ hmax.
% add_point_snap   : snap new vertex to existing grid using tolerance TOL.
% materialize_map_points : turn the snap map back into XY array.
% tri_edges_attach : list unique mesh edges + triangle adjacency (for boundary detection)
% ------------- helper -------------
function A = collect_areas(idx, Lumen, Nc)
    if isempty(idx), A = []; return; end
    A = NaN(numel(idx),1);
    for k = 1:numel(idx)
        ci = idx(k);
        if ci>=1 && ci<=Nc && Lumen(ci).NumRegions>0
            A(k) = area(Lumen(ci));
        end
    end
    A = A(isfinite(A) & A>0);
end

function xy2 = densifyPolyline(xy, hmax)
    xy2 = xy(1,:);
    for i=1:size(xy,1)-1
        p=xy(i,:); q=xy(i+1,:); L=hypot(q(1)-p(1),q(2)-p(2));
        n=max(1, ceil(L/hmax));
        t=linspace(0,1,n+1).';
        seg = p.*(1-t) + q.*t;
        xy2=[xy2; seg(2:end,:)]; %#ok<AGROW>
    end
end

function idx = add_point_snap(p, tol, map)
    px = round(p(1)/tol)*tol;  py = round(p(2)/tol)*tol;
    key = sprintf('%.12g,%.12g', px, py);
    if isKey(map,key), idx = map(key); else, idx = map.Count + 1; map(key) = idx; end
end

function XY = materialize_map_points(map)
    XY = zeros(map.Count,2); kk = map.keys;
    for i=1:numel(kk)
        [px,py] = sscanf(kk{i}, '%f,%f');
        XY(map(kk{i}),:) = px(:).';
    end
end

function [E, E2T] = tri_edges_attach(T)
    allE = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
    allE = sort(allE,2);
    triIx = [(1:size(T,1))'; (1:size(T,1))'; (1:size(T,1))'];
    [E,~,ic] = unique(allE,'rows');
    E2T = zeros(size(E,1),2);
    for k=1:numel(ic)
        e=ic(k); t=triIx(k);
        if E2T(e,1)==0, E2T(e,1)=t; else, E2T(e,2)=t; end
    end
end

%% Save data for further modelling 
save(name_for_saving);
