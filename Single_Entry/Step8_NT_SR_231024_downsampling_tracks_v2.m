%%
%Downsamplng Code - NT(NT_231021_downsampling_tracks.m); Recalculating data - SR

%% The script downsamples all the tracks to the number of ground truth time points
%%% The required data is recalculated and saved

clc
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

ARTfiles = dir(fullfile(sav_path,'*_ART_Track.mat'));
ART1files = dir(fullfile(sav_path,'*_ART_Track1.mat'));
MATfiles = dir(fullfile(sav_path,'*_MAT_16_18_Track1.mat'));
TETfiles = dir(fullfile(sav_path,'*_TET_Track.mat'));

st_inter=150;
end_inter=183;
last_tp=282; %number of phase images before interpolation
last_tp_int = 777; %number of phase images after interpolation 

% ART Tracks
for it0=1:numel(ARTfiles)

    ART_files = fullfile(sav_path, ARTfiles(it0).name);
    load(ART_files);

    endi=size(Mask3,2)-(last_tp-end_inter);
    
    tr1=Mask3(1:st_inter-1);
    tr2=Mask3(st_inter:endi-1);
    tr3=Mask3(endi:end);

ground_truth = cell(0);
numImages = numel(tr2);

for i = 1:numImages-1
    if mod(i, 16) == 1
        ground_truth{end+1} = tr2{i};
    end
end

tr_final = [tr1 ground_truth tr3];
Mask3=tr_final;

ARTname=ARTfiles(it0).name;
ARTparts = strsplit(ARTname, '.');
newARTname = [ARTparts{1} '_DS.mat'];

% Recalculating Data
all_ob = SR_240222_cal_allob(no_obj,Mask3,1:numel(Mask3));

cell_artifacts = [];
cell_exists=zeros(2,size(all_ob,1));   
for itt2 = 1:size(all_ob,1)
    if isempty(find(all_ob(itt2,:)>0,1,'first'))
        cell_artifacts = [cell_artifacts itt2];
    else
        cell_exists(1,itt2)=find(all_ob(itt2,:)>0,1,'first');
        cell_exists(2,itt2)=find(all_ob(itt2,:)>0,1,'last');
    end
end

if ~isempty(cell_artifacts)
    all_ccel = 1:no_obj;
    good_cells = setdiff(all_ccel, cell_artifacts);
    good_cells = sort(good_cells);

    for iv=1:size(good_cells,2)
        for its=1:numel(Mask3)
            pix = find(Mask3{1,its}==good_cells(iv));
            Mask3{1,its}(pix) = iv;
        end
    end
    no_obj = size(good_cells,2);
    all_ob = SR_240222_cal_allob(no_obj,Mask3,1:numel(Mask3));

    cell_exists=zeros(2,size(all_ob,1));   
    for itt2 = 1:size(all_ob,1)
        cell_exists(1,itt2)=find(all_ob(itt2,:)>0,1,'first');
        cell_exists(2,itt2)=find(all_ob(itt2,:)>0,1,'last');
    end
end

save(fullfile(sav_path,newARTname),  'Mask3', 'all_ob', 'ccell2', 'cell_exists', 'no_obj', '-v7.3');

end 

% ART Track1
for it0=1:numel(ART1files)

    ART1_files = fullfile(sav_path, ART1files(it0).name);
    load(ART1_files);

    endi=size(Mask3,2)-(last_tp-end_inter);
    
    tr1=Mask3(1:st_inter-1);
    tr2=Mask3(st_inter:endi-1);
    tr3=Mask3(endi:end);

ground_truth = cell(0);
numImages = numel(tr2);

for i = 1:numImages-1
    if mod(i, 16) == 1
        ground_truth{end+1} = tr2{i};
    end
end

tr_final = [tr1 ground_truth tr3];
Mask3=tr_final;

ART1name=ART1files(it0).name;
ART1parts = strsplit(ART1name, '.');
newART1name = [ART1parts{1} '_DS.mat'];

% Calculating the sizes
all_ob = SR_240222_cal_allob(no_obj,Mask3,1:numel(Mask3));

cell_artifacts = [];
cell_exists=zeros(2,size(all_ob,1));   
for itt2 = 1:size(all_ob,1)
    if isempty(find(all_ob(itt2,:)>0,1,'first'))
        cell_artifacts = [cell_artifacts itt2];
    else
        cell_exists(1,itt2)=find(all_ob(itt2,:)>0,1,'first');
        cell_exists(2,itt2)=find(all_ob(itt2,:)>0,1,'last');
    end
end

if ~isempty(cell_artifacts)
    all_ccel = 1:no_obj;
    good_cells = setdiff(all_ccel, cell_artifacts);
    good_cells = sort(good_cells);

    for iv=1:size(good_cells,2)
        for its=1:numel(Mask3)
            pix = find(Mask3{1,its}==good_cells(iv));
            Mask3{1,its}(pix) = iv;
        end
    end
    no_obj = size(good_cells,2);
    all_ob = SR_240222_cal_allob(no_obj,Mask3,1:numel(Mask3));

    cell_exists=zeros(2,size(all_ob,1));   
    for itt2 = 1:size(all_ob,1)
        cell_exists(1,itt2)=find(all_ob(itt2,:)>0,1,'first');
        cell_exists(2,itt2)=find(all_ob(itt2,:)>0,1,'last');
    end
    
end

save(fullfile(sav_path,newART1name ),  'Mask3', 'all_ob', 'ccell2', 'cell_exists','no_obj', '-v7.3');

end 

% MAT Tracks
for it0=1:numel(MATfiles)

    MAT_files = fullfile(sav_path, MATfiles(it0).name);
    load(MAT_files);

    MATC = Matmasks;
    for itx1 = 1:last_tp_int
        if itx1 > numel(Matmasks)
            MATC{1,itx1} = [];
        end
    end
    Matmasks = MATC;

    endi=size(Matmasks,2)-(last_tp-end_inter);
    
    tr1=Matmasks(1:st_inter-1);
    tr2=Matmasks(st_inter:endi-1);
    tr3=Matmasks(endi:end);

    ground_truth = cell(0);
    numImages = numel(tr2);
    
    for i = 1:numImages-1
        if mod(i, 16) == 1
            ground_truth{end+1} = tr2{i};
        end
    end

    tr_final = [tr1 ground_truth tr3];
    Matmasks=tr_final;
    
    MATname=MATfiles(it0).name;
    MATparts = strsplit(MATname, '.');
    newMATname = [MATparts{1} '_DS.mat'];

% Recalculating Data
all_obj = SR_240222_cal_allob(no_obj,Matmasks,1:numel(Matmasks));

cell_artifacts = [];
cell_data=zeros(size(all_obj,1),2);   
for itt2 = 1:size(all_obj,1)
    if isempty(find(all_obj(itt2,:)>0,1,'first'))
        cell_artifacts = [cell_artifacts itt2];
    else
        cell_data(itt2,1)=find(all_obj(itt2,:)>0,1,'first');
        cell_data(itt2,2)=find(all_obj(itt2,:)>0,1,'last');
    end
end

if ~isempty(cell_artifacts)
    all_ccel = 1:no_obj;
    good_cells = setdiff(all_ccel, cell_artifacts);
    good_cells = sort(good_cells);

    for iv=1:size(good_cells,2)
        for its=1:numel(Matmasks)
            pix = find(Matmasks{1,its}==good_cells(iv));
            Matmasks{1,its}(pix) = iv;
        end
    end
    no_obj = size(good_cells,2);
    all_obj = SR_240222_cal_allob(no_obj,Matmasks,1:numel(Matmasks));

    cell_data=zeros(size(all_obj,1),2);
    for itt2 = 1:size(all_obj,1)
        cell_data(itt2,1)=find(all_obj(itt2,:)>0,1,'first');
        cell_data(itt2,2)=find(all_obj(itt2,:)>0,1,'last');
    end
end

save(fullfile(sav_path,newMATname ), "Matmasks","cell_data","no_obj","shock_period","all_obj" ,'-v7.3');

end 

% TET Tracks
for it0=1:numel(TETfiles)

    TET_files = fullfile(sav_path, TETfiles(it0).name);
    load(TET_files);

    TETC = TETmasks;
    for itx1 = 1:last_tp_int
        if itx1 > numel(TETmasks)
            TETC{1,itx1} = [];
        end
    end
    TETmasks = TETC;

    endi=size(TETmasks,2)-(last_tp-end_inter);
    
    tr1=TETmasks(1:st_inter-1);
    tr2=TETmasks(st_inter:endi-1);
    tr3=TETmasks(endi:end);

ground_truth = cell(0);
numImages = numel(tr2);

for i = 1:numImages-1
    if mod(i, 16) == 1
        ground_truth{end+1} = tr2{i};
    end
end

tr_final = [tr1 ground_truth tr3];
TETmasks=tr_final;

TETname=TETfiles(it0).name;
TETparts = strsplit(TETname, '.');
newTETname = [TETparts{1} '_DS.mat'];

% Recalculating Data
TET_Size = SR_240222_cal_allob(TET_obj,TETmasks,1:numel(TETmasks));

cell_artifacts = [];
TET_exists=zeros(size(TET_Size,1),2);   
for itt2 = 1:size(TET_Size,1)
    if isempty(find(TET_Size(itt2,:)>0,1,'first'))
        cell_artifacts = [cell_artifacts itt2];
    else
        TET_exists(itt2,1)=find(TET_Size(itt2,:)>0,1,'first');
        TET_exists(itt2,2)=find(TET_Size(itt2,:)>0,1,'last');
    end
end

if ~isempty(cell_artifacts)
    all_ccel = 1:TET_obj;
    good_cells = setdiff(all_ccel, cell_artifacts);
    good_cells = sort(good_cells);

    for iv=1:size(good_cells,2)
        for its=1:numel(TETmasks)
            pix = find(TETmasks{1,its}==good_cells(iv));
            TETmasks{1,its}(pix) = iv;
        end
    end
    TET_obj = size(good_cells,2);
    TET_Size = SR_240222_cal_allob(TET_obj,Matmasks,1:numel(TETmasks));

    TET_exists=zeros(size(TET_Size,1),2);   
    for itt2 = 1:size(TET_Size,1)
        TET_exists(itt2,1)=find(TET_Size(itt2,:)>0,1,'first');
        TET_exists(itt2,2)=find(TET_Size(itt2,:)>0,1,'last');
    end
end

save(fullfile(sav_path,newTETname),  'TETmasks', 'shock_period', 'TET_exists', 'tet_masks_exists_tp', 'TET_obj', 'TET_Size', 'thresh', 'thresh_next_cell', 'thresh_perecent', '-v7.3');

end 


