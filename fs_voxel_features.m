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
    
    output.([side(hemisphere) 'h']) = zeros(length(pialv),9);

    TotalArea = NaN(length(pialv_ds),1);
    SmoothArea = NaN(length(pialv_ds),1);
    AvgThickness = NaN(length(pialv_ds),1);
    GaussCurv = NaN(length(pialv_ds),1);
    
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
            
            % Label smooth pial
            label_smooth = label(nearest_opial);            

            % Find the smooth area
            sids = find(label_smooth == 1);
            [liaa] = ismember(opialf,sids);
            aid = (sum(liaa,2) == 3);

            SmoothArea(point) = sum(calcTriangleArea(opialf(aid,:),opialv));

            % Label the ds smooth pial
            label_smooth_ds = label(nearest_opialds);

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
    
    output.lh = array2table(output.lh, 'VariableNames',{'PialArea','SmoothPialArea','PialAreaRaw','SmoothPialAreaRaw', ...
        'AvgCortThickness','GaussCurv','K','I','S'});
    output.rh = array2table(output.rh, 'VariableNames',{'PialArea','SmoothPialArea','PialAreaRaw','SmoothPialAreaRaw', ...
        'AvgCortThickness','GaussCurv','K','I','S'});
end

% Write to file
save([subject '/surf/voxel_features.mat'], 'output')
end