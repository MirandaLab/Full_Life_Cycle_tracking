%% Correcting mating tracks using SpoSeg Tracks

clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';
    
mat_track_path = [sav_path pos '_MAT_16_18_Track'];
file_list = dir([mat_track_path, '*']);
mat = load([sav_path file_list(1).name]); % load mat tracks

if mat.no_obj ~= 0 % number of positive mat detections
                
    shock_period = mat.shock_period;
    MTrack = mat.Matmasks;
    no_obj = mat.no_obj;
    cell_data = mat.cell_data;
    
    art_track_path = [sav_path pos '_ART_Track'];
    file_list = dir([art_track_path, '*']);
    art = load([sav_path file_list(1).name]);
    art_masks  = art.Mask3;
    mat_artifacts = [];
    
    
    for its = 1:numel(MTrack)
        if ~isempty(MTrack{its})
            MTrack{its} = imresize(MTrack{its},[size(art_masks{its})],'nearest');
        end
    end
    
    tp_end = numel(art_masks);
    if numel(MTrack) ~= tp_end
        for its = numel(MTrack)+1:tp_end
            MTrack{its} = uint16(zeros(size(MTrack{min(cell_data(:,1))})));
        end
    end
    % In MatSeg Tracks correction, mating events whose variance of the eccentricity 
    % of each mask across the total valid time points is greater than 0.02 are artifacts
    cor_data = zeros(3,no_obj);
    size_cell = zeros(no_obj,numel(MTrack));
    morph_data = zeros(no_obj,numel(MTrack));
    outlier_tps = cell(1,no_obj);
    good_tps = cell(1,no_obj);
    for iv = 1:no_obj
        int = cell_data(iv,1):cell_data(iv,2);
        for its = int
            M = (MTrack{1,its}==iv);
            size_cell(iv,its) = sum(M(:));
            val = regionprops(M,'Eccentricity');
            morph_data(iv,its) = val(1).Eccentricity;
        end
        cor_data(1,iv) = mean(size_cell(iv,int)); %average
        cor_data(2,iv) = std(size_cell(iv,int)); %standard deviation
        cor_data(3,iv) = 1 * cor_data(2,iv); %threshold
        outlier_tps{iv} = int(abs(size_cell(iv,int)-cor_data(1,iv)) > cor_data(3,iv));
        good_tps{iv} = setdiff(int,outlier_tps{iv});
    end
            
    for iv = 1:no_obj
        int = cell_data(iv,1):cell_data(iv,2);
        if var(morph_data(iv,int)) > 0.02
            mat_artifacts = [mat_artifacts iv];
        end
    end
            
    for iv = 1:no_obj
        outlier = sort(outlier_tps{iv});
        good = sort(good_tps{iv});
        int = cell_data(iv,1):cell_data(iv,2);
        while ~isempty(outlier)
            its = min(outlier);
            if its == min(int)
                gtp = its;
            else
                gtp = max(good(good < its)); %previous good time point
                if isempty(gtp)
                    gtp = min(good(good > its)); %future good time point
                end
            end
            A = uint16(art_masks{1,its}); %figure();imagesc(A)
            M1 = (MTrack{1,gtp}==iv); %figure(1);imagesc(M1)
            M2 = bwmorph(M1,'thin',30); %figure(2);imagesc(M2)
            M3 = A.*uint16(M2); %figure(1);imagesc(M3);pause;
            indx = unique(A(M3~=0))';
            if ~isempty(indx)
                X1 = zeros(size(MTrack{1,its}));
                for itt2 = indx
                    if sum(sum(M3==itt2)) > 5
                        X1(A==itt2) = 1;
                    end
                end
                X1 = imfill(X1, 'holes');
                X2 = bwlabel(logical(X1));
                if max(X2(:))>1
                else
                if abs(sum(X1(:))-cor_data(1,iv)) > 2 * cor_data(2,iv)
                    MTrack{1,its}(MTrack{1,its}==iv) = 0;
                    MTrack{1,its}(MTrack{1,gtp}==iv) = iv;
                else
                    MTrack{1,its}(MTrack{1,its}==iv) = 0;
                    MTrack{1,its}(X1==1) = iv;
                end
                end
            end
            outlier = sort(setdiff(outlier,its));
            good = sort([good its]);
        end
    end
            
    for iv = 1:no_obj
        if cell_data(iv,2) ~= tp_end
            count = 0;
            for its = cell_data(iv,2)+1:tp_end
                A = art_masks{1,its}; %figure(3);imagesc(A)
                M1 = (MTrack{1,its-1}==iv); %figure(1);imagesc(M1)
                M2 = bwmorph(M1,'thin',30); %figure(2);imagesc(M2)
                M3 = uint16(A).*uint16(M2); %figure(1);imagesc(M3);pause;
                indx = unique(A(M3~=0))';
                if ~isempty(indx)
                    X1 = zeros(size(MTrack{1,its}));
                    for itt2 = indx
                        if sum(sum(M3==itt2)) > 5
                            X1(A==itt2) = 1;
                        end
                    end
                    if abs(sum(X1(:))-cor_data(1,iv)) > 2 * cor_data(2,iv)
                        count = count + 1;
                        MTrack{1,its}(MTrack{1,its-1}==iv) = iv;
                    else
                        MTrack{1,its}(X1==1) = iv;
                    end
                else
                    count = count + 1;
                    MTrack{1,its}(MTrack{1,its-1}==iv) = iv;
                end
            end
            if count/numel(cell_data(iv,1):tp_end) > 0.8
                mat_artifacts = [mat_artifacts iv];
            end
        end
    end

    % Remove cell artifacts and rename
    if ~isempty(mat_artifacts)
        all_ccel = 1:no_obj;
        mat_artifacts = unique(mat_artifacts);
        for iv = mat_artifacts
            for its = 1:numel(MTrack)
                pix = find(MTrack{1,its}==iv);
                MTrack{1,its}(pix) = 0;
            end
        end
            
        good_cells = setdiff(all_ccel, mat_artifacts);
        good_cells = sort(good_cells);
        
        for iv=1:size(good_cells,2)
            for its=1:numel(MTrack)
                pix = find(MTrack{1,its}==good_cells(iv));
                MTrack{1,its}(pix) = iv;
            end
        end
        no_obj = numel(good_cells);
    end
        
    % Recalculating MAT Data
    all_obj_new = SR_240222_cal_allob(no_obj,MTrack,1:numel(MTrack));
    cell_data_new = SR_240222_cal_celldata(all_obj_new,no_obj);

    cell_data = cell_data_new;
    all_obj = all_obj_new;
    Matmasks = MTrack;

    save([sav_path pos '_MAT_16_18_Track1.mat'],"Matmasks","all_obj","cell_data","no_obj","shock_period","mat_artifacts","-v7.3");

end