%% THIS IS FOR MAT TRACK 1 TRACKS

clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

tet_track_path = [sav_path pos '_TET_Track_DS'];
file_list = dir([tet_track_path, '*']);
load([sav_path file_list(1).name]);  % load tet track

tet_track_path = [sav_path pos '_TET_ID'];
file_list = dir([tet_track_path, '*']);
load([sav_path file_list(1).name]); % load tet IDs

desc_path = [sav_path pos '_descendants_new_art'];
file_list = dir([desc_path, '*']);
load([sav_path file_list(1).name]); % load descendants information

mat_track_path = [sav_path pos '_MAT_16_18_Track1_DS'];
file_list = dir([mat_track_path, '*']);
load([sav_path file_list(1).name]); % Load mat track

if no_obj ~= 0 % positive number of MAT Detections
    masks = TETmasks;
    MTrack = Matmasks;
    shock_period = shock_period;
    for i = 1:numel(masks)
        if ~isempty(masks{i})
        masks{i} = imresize(masks{i},size(I3),"nearest");
        end
    end
    
    start = shock_period(1,2)+1;
    tp_end = numel(MTrack); %190;
    % nutr_stress_tp = numel(MTrack);%190;
    % intx2 = start:nutr_stress_tp;
    int = min(cell_data(:,1)):max(cell_data(:,2)); %cell_data(:,1)
    
    % Checking which cell from mat tracks belongs to which region created using watershed from the previous step

    region = [];
    amt = []; 
    k = 0;
    for iv = 1:no_obj
        I12 = zeros(size(uint16(MTrack{1,min(int)})));
        kx = 0;
        for its = int
            I11 = logical(MTrack{1,its}==iv);
            if sum(I11(:)) > 0
                kx = kx+1;
                if kx >= 1 && kx <= 2 %kx==1 % the first 2 locations of the cells are given more weightage
                    I11 = 1000.*I11;%weight.*I11; %figure;imagesc(I11)
                end
            end
            I12 = I12 + I11;
        end
        I13 = uint16(logical(I12)).*uint16(I3);
        pix = unique(I13)';
        pix = pix(pix~=0);
        if ~isempty(pix)
            for i = 1:size(pix,2)
                amt(iv,i) = sum(sum(I13(:)==pix(i)));
            end
            k = k + 1;
            [val,ind] = max(amt(iv,:));
            region(k,1) = iv; % cell number in mat tracks
            region(k,2) = pix(ind); %region
        end
    end
    
    unique_regions = unique(region(:, 2)); %% All cells in artilife that belong to a particular sector
    cell_arrays = cell(length(unique_regions), 1);
    for i = 1:length(unique_regions)
        current_region = unique_regions(i);
        cells_in_current_region = region(region(:, 2) == current_region,1);
        cell_arrays{unique_regions(i)} = uint16(cells_in_current_region); % 
    end
    
    if ~isempty(alive_tets)
        cell_re = {};
        TET_ind = zeros(numel(unique_regions),3);
        common_indices = cell(1,TET_obj);
        amt_1 = [];
        for iv = 1:TET_obj % loop through tet indicies
            if TET_ID(1,iv) ~= -1
                if ~isempty(find(alive_tets == iv,1))
            if TET_exists(iv,2) >= shock_period(1,2)+1
                T1 = masks{1,shock_period(1,2)+1}==iv;
            else
                T1 = masks{1,TET_exists(iv,2)}==iv;
            end
            %T1 = masks{1,TET_exists(iv,1)}==iv; % figure;imagesc(T1);
            T2 = uint16(I3).*uint16(T1);
            TET_ind(iv,1) = iv; % TET number
            TET_ind(iv,3) = TET_ID(1,iv);%max(T3(:)); % TET index in artilife %TET_ID !! USE TET ID INSTEAD
    
            pix = unique(T2)';
            pix = pix(pix~=0);
            if ~isempty(pix)
            for i = 1:size(pix,2)
                amt_1(iv,i) = sum(sum(T2(:)==pix(i)));
            end
            [val,ind] = max(amt_1(iv,:));
            TET_ind(iv,2) = pix(ind); % TET region
            end
                
                else
                    TET_ind(iv,1) = iv; % TET number
                    TET_ind(iv,3) = TET_ID(1,iv);
                end
            end
        end
    
        for ixx = 1:size(cell_arrays,1)
            cell_arrays{ixx,1}
            tet_no = find(TET_ind(:,2)==ixx);
            descendants_data{tet_no,4} = cell_arrays{ixx,1};
        end
    end
    descendants_data
    save([sav_path pos '_final_descendants.mat'],"I3","descendants_data","alive_tets","TET_obj","-v7.3");
end
