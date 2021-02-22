function [output] = fs_voxel_features(subject)

% Subject could be e.g. '/nobackup/b9014486/data/LesionProject/D021'
% Needs output table from extract_FreeSurferHemi_features.m function

libdir = '../FreesurferExtractFeatures/lib';
addpath(libdir)
addpath([libdir '/FSmatlab/'])
addpath([libdir '/iso2mesh/'])


side = 'lr';

output = struct;

% Order of variables in output
% At_dash (using GC CH)
% Ae_dash (using GC CH)
% At raw
% Ae raw
% AvgThickness
% GaussCurv
% GaussCurvPial
% GaussCurvSmoothPial
% K
% I
% S

% Spherical radius around point (mm)
% r = 15;


for hemisphere = 1:2
% for hemisphere = hs
    
    if hemisphere == 1
        hemi = "left";
    else
        hemi = "right";
    end
    
    % Subject ID
    [~,name] = fileparts(subject);
    
    % Get hemi k (not logged), table is output from extract_FreeSurferHemi_features.m
%     tbl = load('C:/linux/freesurfer_subjects/LesionProject/LesionProject_hemi.mat')
    tbl = load('../LesionProject_hemi.mat');
    tbl = tbl.tbl;
    k_hemi = tbl((tbl.SubjectID == string(name)) & (tbl.Hemisphere == hemi),:).PialArea * ...
        sqrt(tbl((tbl.SubjectID == string(name)) & (tbl.Hemisphere == hemi),:).AvgCortThickness) * ...
        tbl((tbl.SubjectID == string(name)) & (tbl.Hemisphere == hemi),:).SmoothPialArea^(-5/4);

    % Path to subject's files
    pathpre = [subject '/surf/' side(hemisphere)];

    [thickness, ~]  = read_curv([pathpre, 'h.thickness']);
    [pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial']);
    [opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed']);
%     [spherev,spheref] = freesurfer_read_surf([pathpre, 'h.sphere']);
    
    
    output.([side(hemisphere) 'h']) = zeros(length(pialv),11);

    TotalArea = NaN(length(pialv_ds),1);
    SmoothArea = NaN(length(pialv_ds),1);
    AvgThickness = NaN(length(pialv_ds),1);
    GaussCurv = NaN(length(pialv_ds),1);
    GaussCurvPial = NaN(length(pialv_ds),1);
    GaussCurvSmoothPial = NaN(length(pialv_ds),1);
%     GaussCurvSphere= NaN(length(pialv_ds),1);
    
    % Downsample pial surface
    [pialv_ds,pialf_ds] = meshresample(pialv,pialf,0.05);

    % Downsample smooth pial
    [opialv_ds, opialf_ds] = meshresample(opialv,opialf,0.1);
    
    % Find nearest point for each point in pial (and at the same time also sphere)
    % Split into two blocks so the size of the matrix doesn't crash the
    % memory
    pialv_ds = single(pialv_ds);
    d1 = pialv_ds(:,1) - pialv(1:floor(size(pialv,1)/2),1)';
    d2 = pialv_ds(:,2) - pialv(1:floor(size(pialv,1)/2),2)';
    d3 = pialv_ds(:,3) - pialv(1:floor(size(pialv,1)/2),3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_sphere] = min(d);
    clear d
    
    d1 = pialv_ds(:,1) - pialv((floor(size(pialv,1)/2)+1):end,1)';
    d2 = pialv_ds(:,2) - pialv((floor(size(pialv,1)/2)+1):end,2)';
    d3 = pialv_ds(:,3) - pialv((floor(size(pialv,1)/2)+1):end,3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_sphere2] = min(d);
    clear d
    nearest_sphere = [nearest_sphere nearest_sphere2];
    clear nearest_sphere2
    
    
    % Other way around (for each point on ds pial the point in pial that is
    % nearest)
    d1 = pialv(:,1) - pialv_ds(1:floor(size(pialv_ds,1)/2),1)';
    d2 = pialv(:,2) - pialv_ds(1:floor(size(pialv_ds,1)/2),2)';
    d3 = pialv(:,3) - pialv_ds(1:floor(size(pialv_ds,1)/2),3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_pialds] = min(d);
    clear d
    
    d1 = pialv(:,1) - pialv_ds((floor(size(pialv_ds,1)/2)+1):end,1)';
    d2 = pialv(:,2) - pialv_ds((floor(size(pialv_ds,1)/2)+1):end,2)';
    d3 = pialv(:,3) - pialv_ds((floor(size(pialv_ds,1)/2)+1):end,3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_pialds2] = min(d);
    clear d
    nearest_pialds= [nearest_pialds nearest_pialds2];
    clear nearest_pialds2
    
    % Nearest point for points in smooth pial
    pialv_ds = single(pialv_ds);
    d1 = pialv_ds(:,1) - opialv(:,1)';
    d2 = pialv_ds(:,2) - opialv(:,2)';
    d3 = pialv_ds(:,3) - opialv(:,3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opial] = min(d);

    % Nearest point for points in ds smooth pial
    d1 = pialv_ds(:,1) - opialv_ds(:,1)';
    d2 = pialv_ds(:,2) - opialv_ds(:,2)';
    d3 = pialv_ds(:,3) - opialv_ds(:,3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opialds] = min(d);
    clear d
    
    % Gauss curv on ds pial
%     GC = getGaussianCurvPart(pialf_ds, pialv_ds, 1:length(pialv_ds));
    
    % Gauss curv on ds smooth pial
%     GCsmooth = getGaussianCurvPart(opialf_ds, opialv_ds, 1:length(opialv_ds));
    
    % Gauss curv on sphere
%     GCsphere = getGaussianCurvPart(spheref, spherev, 1:length(spherev));

    % Find thickness of vertices of pial closest to each ds vertex
    thickness_ds = thickness(nearest_pialds);

    % Thickness of faces in ds pial
    ThicknessFB = makeFacebased(thickness_ds, pialf_ds);

    
    for point = 1:length(pialv_ds)

        % Check if it is on the midline or if the thickness is too small
        % (for healthy adults)
        if thickness_ds(point) < 1.8 || abs(pialv_ds(point,1)) < 10
            continue;
        else

            % Index of point being looked at
            neighbourIDs = point;

            % Coordinates of that point
            vertex = pialv_ds(neighbourIDs,:);

            % Radius needed for patch dependent on thickness at that point
            r = thickness_ds(point)/(sqrt(pi)*k_hemi^2); % ~ 11-30mm
%             r = thickness_ds(point)/(k_hemi^2);

            % Repeat until new neighbours are all > r mm from the point
            rmin = 0;
            while rmin < r

                % Find neighbours of neighbours
                fid = ismember(pialf_ds(:,1),neighbourIDs) | ...
                    ismember(pialf_ds(:,2),neighbourIDs) | ...
                    ismember(pialf_ds(:,3),neighbourIDs);

                % Check how far all new ones are from vertex, add the ones that are
                % close enough
                new = setdiff(unique(pialf_ds(fid,:)), neighbourIDs);
                d = sqrt(sum((pialv_ds(new,:)-vertex).^2,2));
                rmin = min(d);

                neighbourIDs = [neighbourIDs; new(d<r)];
            end


            % Indices of faces that contain neighbour vertices
%             fid = find(ismember(pialf_ds(:,1),neighbourIDs) & ...
%                 ismember(pialf_ds(:,2),neighbourIDs) & ...
%                 ismember(pialf_ds(:,3),neighbourIDs));

            % Total area
            % Only use area of faces that are completely within circle
            [liaa] = ismember(pialf_ds, neighbourIDs);
            aid = (sum(liaa,2) == 3);

            TotalAreai = zeros(length(aid), 1);
            TotalAreai(aid) = calcTriangleArea(pialf_ds(aid,:), pialv_ds);
            TotalArea(point) = sum(TotalAreai);

            % Average thickness in circle
            AvgThickness(point) = sum(ThicknessFB(aid>0).*TotalAreai(aid>0))/TotalArea(point);

            % Label ds pial
            label = zeros(length(pialv_ds),1);
            label(neighbourIDs) = 1;

            % Label pial/sphere
%             label_sphere = label(nearest_sphere);
            
            % Get curvature for patch on sphere
%             GaussCurvSphere(point) = sum(GCsphere(find(label_sphere)));
            
            % Label smooth pial
%             label_smooth = matchSurfLabel(label, pialv_ds, opialv);
            label_smooth = label(nearest_opial);            

            % Find the smooth area
            sids = find(label_smooth == 1);
            [liaa] = ismember(opialf,sids);
            aid = (sum(liaa,2) == 3);

            SmoothArea(point) = sum(calcTriangleArea(opialf(aid,:),opialv));

            % Label the ds smooth pial
%             label_smooth_ds = matchSurfLabel(label,pialv_ds,opialv_ds);
            label_smooth_ds = label(nearest_opialds);

%             GaussCurvPial(point) = sum(GC(neighbourIDs));

            ov_ids = find(label_smooth_ds == 1);

%             GaussCurvSmoothPial(point) = sum(GCsmooth(ov_ids));

            Lobepoints = opialv_ds(ov_ids,:);

            if length(Lobepoints) > 4

                CHSf_ds = convhull(Lobepoints);

                Gru = getGaussianCurvPart(CHSf_ds,Lobepoints,1:length(ov_ids));

                isBoundary = zeros(length(ov_ids),1);

                for kl = 1:length(ov_ids)
                    fid = opialf_ds(:,1)==ov_ids(kl) | ...
                          opialf_ds(:,2)==ov_ids(kl) | ...
                          opialf_ds(:,3)==ov_ids(kl);

                    v_neigh = unique(opialf_ds(fid,:));
                    ncolors = numel(unique(label_smooth_ds(v_neigh)));
                    if ncolors > 1
                        isBoundary(kl) = 1;
                    end
                end

                GaussCurv(point)  = sum(Gru(isBoundary == 0));
            end
        end
    end
    
    % Correct surface areas
    At_dash = TotalArea.*4*pi./GaussCurv;
    Ae_dash = SmoothArea.*4*pi./GaussCurv;
    
    % Determine cutoff point for CH curvature (Affects At_dash, Ae_dash)
    slope = NaN(1,1000);
    cutoff = 0.001:0.001:1;
    
    % Find point at which the slope stabilises
    for k = 1:1000
        
        % Determine which points to use
        usepoints1 = find(GaussCurv >= cutoff(k));

        % Correct At and Ae, then get X and Y
        At_dash1 = TotalArea(usepoints1).*(4*pi)./GaussCurv(usepoints1);
        Ae_dash1 = SmoothArea(usepoints1).*(4*pi)./GaussCurv(usepoints1);

        X1 = log10(abs(Ae_dash1));
        Y1 = log10(abs(At_dash1).*sqrt(AvgThickness(usepoints1)));

        if sum(X1 ~= Inf) > 1
            m = fitlm(X1(X1 ~= Inf),Y1(X1 ~= Inf));
            slope(1,k) = m.Coefficients.Estimate(2);
        end
    end
    
    change = zeros(1,900);
    % Find cutoff after which the slope stabilises
    for k = 1:100
        change = change + abs(slope(1,1:900) - slope(1,(1+k):(900+k)));
    end
    
    change = change./100;
    change = change < 0.002*slope(1,1:900);
    
    stable = find(change, 1, 'first');
    
    if ~isempty(stable)
        % Set values where CH GC is below stable point and points on midline to NaN
        At_dash(GaussCurv < stable*0.001) = NaN;
        Ae_dash(GaussCurv < stable*0.001) = NaN;
    else
        disp(strcat("Slope did not stabilise in ", string(side(hemisphere)), "h."))
        disp(" ")
    end
    
    K = log10(At_dash) - 5/4*log10(Ae_dash) + ...
        1/2*log10(AvgThickness);
    I = log10(At_dash) + log10(Ae_dash) + ...
        2*log10(AvgThickness);
    S = 3/2*log10(At_dash) + 3/4*log10(Ae_dash) - ...
        2*9/4*log10(AvgThickness); 

    % Go back to pial surface
    
    % At_dash (using GC CH)
    % Ae_dash (using GC CH)
    % At raw
    % Ae raw
    % AvgThickness
    % GaussCurv
    % GaussCurvPial
    % GaussCurvSmoothPial
    % K
    % I
    % S
    output.([side(hemisphere) 'h'])(:,1) = At_dash(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,2) = Ae_dash(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,3) = TotalArea(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,4) = SmoothArea(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,5) = AvgThickness(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,6) = GaussCurv(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,7) = GaussCurvPial(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,8) = GaussCurvSmoothPial(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,9) = K(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,10) = I(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,11) = S(nearest_sphere);
    
    output.([side(hemisphere) 'h'])(output.([side(hemisphere) 'h']) == Inf | output.([side(hemisphere) 'h']) == -Inf) = NaN;
    
    
    % Write to files
    fnum = length(pialf);
    write_curv([pathpre 'h.PialArea'], output.([side(hemisphere) 'h'])(:,1), fnum);
    write_curv([pathpre 'h.SmoothPialArea'], output.([side(hemisphere) 'h'])(:,2), fnum);
    write_curv([pathpre 'h.PialAreaRaw'], output.([side(hemisphere) 'h'])(:,3), fnum);
    write_curv([pathpre 'h.SmoothPialAreaRaw'], output.([side(hemisphere) 'h'])(:,4), fnum);
    write_curv([pathpre 'h.AvgCortThickness'], output.([side(hemisphere) 'h'])(:,5), fnum);
    write_curv([pathpre 'h.GaussCurv'], output.([side(hemisphere) 'h'])(:,6), fnum);
    write_curv([pathpre 'h.GaussCurvPial'], output.([side(hemisphere) 'h'])(:,7), fnum);
    write_curv([pathpre 'h.GaussCurvSmoothPial'], output.([side(hemisphere) 'h'])(:,8), fnum);
    write_curv([pathpre 'h.K'], output.([side(hemisphere) 'h'])(:,9), fnum);
    write_curv([pathpre 'h.I'], output.([side(hemisphere) 'h'])(:,10), fnum);
    write_curv([pathpre 'h.S'], output.([side(hemisphere) 'h'])(:,11), fnum);
end
end

