function [output] = fs_voxel_features(subject, libdir, iso2meshdir)
%
% Input:
% - subject: char, path to FreeSurfer folder
% - libdir: char, path to lib-folder: https://github.com/cnnp-lab/CorticalFoldingAnalysisTools/tree/master/lib
% - iso2meshdir: char, path to ISO2MESH libaray: http://iso2mesh.sourceforge.net/
%   Licence: CC-BY
%
% Output: table for each hemisphere, contains values at each point on the
% pial surface
% - At (pial area) corrected using Gaussian curvature of convex hull
% - Ae (smooth pial area) corrected using Gaussian curvature of convex hull
% - At (pial area) raw
% - Ae (smooth pial area) raw
% - Average cortical thickness
% - Gaussian curvature
% - K
% - I
% - S

addpath(libdir)
addpath([libdir '/FSmatlab/'])
addpath(iso2meshdir)

side = 'lr';

output = struct;

% Spherical radius around point (mm)
r = 25;

for hemisphere = 1:2

    % Path to subject's files
    pathpre = [subject '/surf/' side(hemisphere)];

    [thickness, ~]  = read_curv([pathpre, 'h.thickness']);
    [pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial']);
    [opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed']);
    [~,labelDK,~] = read_annotation([subject, '/label/', side(hemisphere), 'h.aparc.annot']);
    
    output.([side(hemisphere) 'h']) = zeros(length(pialv),9);
    
    % Downsample pial surface
    [pialv_ds,pialf_ds] = meshresample(pialv,pialf,0.05);

    % Downsample smooth pial
    [opialv_ds, opialf_ds] = meshresample(opialv,opialf,0.1);

    clear opialv opialf

    TotalArea = NaN(length(pialv_ds),1);
    SmoothArea = NaN(length(pialv_ds),1);
    AvgThickness = NaN(length(pialv_ds),1);
    GaussCurv = NaN(length(pialv_ds),1);
    
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

    % Only keep points on ds pial which are ot on the CC
    label_ds = labelDK(nearest_pialds);
    pialf_ds = pialf_ds((label_ds(pialf_ds(:,1)) ~= 0) & (label_ds(pialf_ds(:,2)) ~= 0) & ...
        (label_ds(pialf_ds(:,3)) ~= 0),:);

    % Points over which to iterate
    usepoints = unique(pialf_ds);
    
    % Label the opial_ds with DK
    % for each point on ds opial the point in pial that is
    % nearest
    d1 = pialv(:,1) - opialv_ds(1:floor(size(opialv_ds,1)/4),1)';
    d2 = pialv(:,2) - opialv_ds(1:floor(size(opialv_ds,1)/4),2)';
    d3 = pialv(:,3) - opialv_ds(1:floor(size(opialv_ds,1)/4),3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opialds] = min(d);
    clear d
    
    d1 = pialv(:,1) - opialv_ds((floor(size(opialv_ds,1)/4)+1):(floor(size(opialv_ds,1)*2/4)),1)';
    d2 = pialv(:,2) - opialv_ds((floor(size(opialv_ds,1)/4)+1):(floor(size(opialv_ds,1)*2/4)),2)';
    d3 = pialv(:,3) - opialv_ds((floor(size(opialv_ds,1)/4)+1):(floor(size(opialv_ds,1)*2/4)),3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opialds2] = min(d);
    clear d
    nearest_opialds= [nearest_opialds nearest_opialds2];
    clear nearest_opialds2
    
    d1 = pialv(:,1) - opialv_ds((floor(size(opialv_ds,1)*2/4)+1):(floor(size(opialv_ds,1)*3/4)),1)';
    d2 = pialv(:,2) - opialv_ds((floor(size(opialv_ds,1)*2/4)+1):(floor(size(opialv_ds,1)*3/4)),2)';
    d3 = pialv(:,3) - opialv_ds((floor(size(opialv_ds,1)*2/4)+1):(floor(size(opialv_ds,1)*3/4)),3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opialds2] = min(d);
    clear d
    nearest_opialds= [nearest_opialds nearest_opialds2];
    clear nearest_opialds2
    
    d1 = pialv(:,1) - opialv_ds((floor(size(opialv_ds,1)*3/4)+1):end,1)';
    d2 = pialv(:,2) - opialv_ds((floor(size(opialv_ds,1)*3/4)+1):end,2)';
    d3 = pialv(:,3) - opialv_ds((floor(size(opialv_ds,1)*3/4)+1):end,3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opialds2] = min(d);
    clear d
    nearest_opialds= [nearest_opialds nearest_opialds2];
    clear nearest_opialds2

    % Only keep point on ds opial which are not on the CC
    label_ods = labelDK(nearest_opialds);
    opialf_ds = opialf_ds((label_ods(opialf_ds(:,1)) ~= 0) & (label_ods(opialf_ds(:,2)) ~= 0) & ...
        (label_ods(opialf_ds(:,3)) ~= 0),:);

    % Nearest point for points in ds smooth pial
    d1 = pialv_ds(:,1) - opialv_ds(:,1)';
    d2 = pialv_ds(:,2) - opialv_ds(:,2)';
    d3 = pialv_ds(:,3) - opialv_ds(:,3)';
    d = d1.^2 + d2.^2 + d3.^2;
    clear d1 d2 d3
    d = sqrt(d);
    [~,nearest_opialds] = min(d);
    clear d

    % Find thickness of vertices of pial closest to each ds vertex
    thickness_ds = thickness(nearest_pialds);
    thickness_ds_smoothed = NaN(length(pialv_ds),1);

    % Smoothing thickness across verices (only used to exclude points)
    for index = 1:length(usepoints)
        point = usepoints(index);
        neighbourIDs = point;
        vertex = pialv_ds(neighbourIDs,:);
        rmin = 0;
        while rmin < 25
            % Find neighbours of neighbours
            fid = ismember(pialf_ds(:,1),neighbourIDs) | ...
                ismember(pialf_ds(:,2),neighbourIDs) | ...
                ismember(pialf_ds(:,3),neighbourIDs);

            new = setdiff(unique(pialf_ds(fid,:)), neighbourIDs);
            d = sqrt(sum((pialv_ds(new,:)-vertex).^2,2));
            rmin = min(d);
            new = new(:);
            neighbourIDs = [neighbourIDs; new];
        end
        thickness_ds_smoothed(point) = mean(thickness_ds(neighbourIDs), 'omitnan');
    end

    % Thickness of faces in ds pial
    ThicknessFB = makeFacebased(thickness_ds, pialf_ds);
    
    for index = 1:length(usepoints)

        point = usepoints(index);

        % Check if it is on the midline or if the thickness is too small
        % (for healthy adults)
        if thickness_ds_smoothed(point) < 1.8 || label_ds(point) == 0
            continue;
        else

            % Index of point being looked at
            neighbourIDs = point;

            % Coordinates of that point
            vertex = pialv_ds(neighbourIDs,:);

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
                new = new(:);

                neighbourIDs = [neighbourIDs; new(d<r)];
            end
            
            % Label ds pial
            label = zeros(length(pialv_ds),1);
            label(neighbourIDs) = 1;
            
            % Find edges of patch
            isBoundary = zeros(length(neighbourIDs),1);
            for kl = 1:length(neighbourIDs)
                fid = pialf_ds(:,1)==neighbourIDs(kl) | ...
                      pialf_ds(:,2)==neighbourIDs(kl) | ...
                      pialf_ds(:,3)==neighbourIDs(kl);

                v_neigh = unique(pialf_ds(fid,:));
                ncolors = numel(unique(label(v_neigh)));
                if ncolors > 1
                    isBoundary(kl) = 1;
                end
            end
            edge = neighbourIDs(isBoundary==1);
            
            % Fill potential holes            
            if ~isempty(edge)
                newPoints = edge(1);
                points_left = 1;
                while points_left
                    % Check we don't start with a vertex that is between a hole and the
                    % edge
                    if ismember(newPoints(1),edge)
                        bndNeigh =  ismember(pialf_ds(:,1),newPoints(1)) | ...
                            ismember(pialf_ds(:,2),newPoints(1)) | ...
                            ismember(pialf_ds(:,3),newPoints(1));

                        bndNeigh = unique(pialf_ds(bndNeigh,:));
                        bndNeigh = bndNeigh(ismember(bndNeigh, neighbourIDs(isBoundary == 1)));
                        if length(bndNeigh) > 3
                            edge(1) = [];
                            if isempty(edge)
                               points_left = 0;
                            else
                               newPoints = edge(1);
                            end
                            continue
                        end
                    end

                    % Find neighbours of neighbours
                    fid = ismember(pialf_ds(:,1),newPoints) | ...
                        ismember(pialf_ds(:,2),newPoints) | ...
                        ismember(pialf_ds(:,3),newPoints);

                    new = unique(pialf_ds(fid,:));
                    edge(ismember(edge,new)) = [];
                    new = setdiff(new, neighbourIDs);
                    new = setdiff(new, newPoints);
                    new = new(:);
                    newPoints = [newPoints; new];
                    newPoints = setdiff(newPoints, neighbourIDs);
                    d = sqrt(sum((pialv_ds(newPoints,:)-vertex).^2,2));

                    if isempty(new)
                        neighbourIDs = [neighbourIDs; newPoints];
                        isBoundary = [isBoundary; zeros(length(newPoints),1)];
                        if isempty(edge)
                           points_left = 0;
                        else
                           newPoints = edge(1);
                        end
                    end

                    if max(d) > 50
                        % Run only on edge until no new points
                        newEdge = 1;
                        while ~isempty(newEdge)
                            fid = ismember(pialf_ds(:,1),newPoints) | ...
                                ismember(pialf_ds(:,2),newPoints) | ...
                                ismember(pialf_ds(:,3),newPoints);

                            new = unique(pialf_ds(fid,:));

                            newEdge = new(ismember(new, edge));
                            edge(ismember(edge,new)) = [];

                            new = setdiff(new, neighbourIDs);
                            new = setdiff(new, newPoints);
                            new = new(:);
                            newPoints = [newPoints; new];
                        end
                        if isempty(edge)
                           points_left = 0;
                        else
                           newPoints = edge(1);
                        end
                    end
                end
            end

            % Relabel ds pial
            label = zeros(length(pialv_ds),1);
            label(neighbourIDs) = 1;

            % Total area
            % Only use area of faces that are completely within circle
            [liaa] = ismember(pialf_ds, neighbourIDs);
            aid = (sum(liaa,2) == 3);

            TotalAreai = zeros(length(aid), 1);
            TotalAreai(aid) = calcTriangleArea(pialf_ds(aid,:), pialv_ds);
            TotalArea(point) = sum(TotalAreai);

            % Average thickness in circle
            AvgThickness(point) = sum(ThicknessFB(aid>0).*TotalAreai(aid>0))/TotalArea(point);
            
%             % Label smooth pial
%             label_smooth = label(nearest_opial);            
% 
%             % Find the smooth area
%             sids = find(label_smooth == 1);
%             [liaa] = ismember(opialf,sids);
%             aid = (sum(liaa,2) == 3);
% 
%             SmoothArea(point) = sum(calcTriangleArea(opialf(aid,:),opialv));

            % Label the ds smooth pial
            label_smooth_ds = label(nearest_opialds);
            
            % Find the smooth area
            sids = find(label_smooth_ds == 1);
            [liaa] = ismember(opialf_ds,sids);
            aid = (sum(liaa,2) == 3);

            SmoothArea(point) = sum(calcTriangleArea(opialf_ds(aid,:),opialv_ds));

            % Gauss. curv. of patch
            ov_ids = find(label_smooth_ds == 1);
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
    
    % Set values where CH GC is below 0.16 to NaN
    At_dash(GaussCurv < 0.16) = NaN;
    Ae_dash(GaussCurv < 0.16) = NaN;
    
    K = log10(At_dash) - 5/4*log10(Ae_dash) + ...
        1/2*log10(AvgThickness);
    I = log10(At_dash) + log10(Ae_dash) + ...
        2*log10(AvgThickness);
    S = 3/2*log10(At_dash) + 3/4*log10(Ae_dash) - ...
        2*9/4*log10(AvgThickness);

    % Go back to pial surface
    output.([side(hemisphere) 'h'])(:,1) = At_dash(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,2) = Ae_dash(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,3) = TotalArea(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,4) = SmoothArea(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,5) = AvgThickness(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,6) = GaussCurv(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,7) = K(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,8) = I(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,9) = S(nearest_sphere);
    
    output.([side(hemisphere) 'h'])(output.([side(hemisphere) 'h']) == Inf | output.([side(hemisphere) 'h']) == -Inf) = NaN;
    
    output.([side(hemisphere) 'h']) = array2table(output.([side(hemisphere) 'h']), 'VariableNames',{'PialArea','SmoothPialArea','PialAreaRaw','SmoothPialAreaRaw', ...
        'AvgCortThickness','GaussCurv','K','I','S'});
end

% Write to file
save([subject '/surf/voxel_features.mat'], 'output')
end