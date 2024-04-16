
clc;
clear;
close all;

thresh_perecent = 0.015;% if number of appearances is less than thresh_percent of total vaild timepoints, the cell is removed
thresh_remove_last_mask = 10; % Number of time points to distinguish germination point from segmentation artifacts where cells are detected as tetrads due to close proximity
thresh_next_cell = 400; % Number of time points to distinguish between same/close enough co-ordinates
thresh = 80; % percentage of times the cell disappears - used in removing artifacts 
shock_period = [122,134]; % time points that need to be removed from tracking

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\']; % Path to segment SpoSeg masks
sav_path = 'E:\SR_Tracking\toy_data\Tracks\'; % Path to save Tracks

path_dir=dir(fullfile(path,'*_Ph3_000_TET_masks.tif'));
fileNames = {path_dir.name};
fileNumbers = zeros(size(fileNames));

for i = 1:numel(fileNames)
    [~, name, ~] = fileparts(fileNames{i});
    startIndex = strfind(name, 'img_') + length('img_');
    numberStr = name(startIndex:end);
    fileNumbers(i) = sscanf(numberStr, '%d');
end

[sortedNumbers, sortedIndices] = sort(fileNumbers);

tet_masks_path = {};
for i = sortedIndices
    img_path = fullfile(path, path_dir(i).name);
    tet_masks_path = [tet_masks_path, img_path];
end

tet_masks = {};
for i = 1:numel(fileNames)
    % +1 because the images start from 0
    tet_masks{sortedNumbers(i)+1} = imread(tet_masks_path{i});
end

% Save all the valid SpoSeg masks in a variable called tet_masks
for i = min(sortedNumbers)+1:size(tet_masks,2)
    if isempty(tet_masks{i})
        tet_masks{i} = uint16(zeros(size(tet_masks{min(sortedNumbers)+1})));
    end
end

%% Remove shock induced timepoints

if ~isempty(shock_period)
    for j = 1:size(shock_period,1)
        for i = shock_period(j,1):shock_period(j,2)
            tet_masks{i} = [];    
        end
    end
end

start = 0;
for its=1:size(tet_masks,2)
    M = uint16(tet_masks{1,its});
    if sum(sum(M)) > 0
        start = its;
        break;
    end
end

%% Tracking all the detections

if start ~=0
    rang = start:numel(tet_masks);
    I2 = tet_masks{start};
    A = zeros(size(tet_masks{start}));
else
    rang = 1:numel(tet_masks);
    I2 = tet_masks{1};
    A = zeros(size(tet_masks{1}));
end

IS6=zeros(size(I2));
TETC = cell(2,numel(tet_masks));
xx = start;
rang2 = rang;
ccel = 1;

while xx ~= 0
    k = 0;
    for im_no=rang2
        if ccel == 1
            I2 = tet_masks{im_no};
        else
            I2 = TETC{2,im_no};
        end
        if isempty(I2)
            continue;
        else
            if im_no==min(rang2)
                ind1 = unique(I2);
                ind1 = ind1(2:end);
                I3=(I2==ind1(1));
                I3A = I3;
            else
                I3A = IS6;
            end
            I3A = bwskel(logical(I3A));

            I2A = I2;

            I3B =uint16(I3A).*uint16(I2A);
            ind = mode(I3B(I3B~=0));
            if ind==0 && ccel == 1
                k = k+1;
                if k > thresh_next_cell
                    for im_no_1 = im_no:rang(end)
                        if ~isempty(tet_masks{im_no_1})
                            TETC{1,im_no_1}=zeros(size(tet_masks{start}));
                        end
                        TETC{2,im_no_1}=tet_masks{im_no_1};
                    end
                    break;
                else
                    TETC{1,im_no}=I3B;
                    TETC{2,im_no}=I2A;
                    continue;
                end
            elseif ind==0 && ccel~=1
                k = k+1;
                if k > thresh_next_cell
                    break;
                else
                    continue;
                end
            end
            k = 0;
            pix=find(I2A==ind);
            pix0=find(I2A~=ind);

            I2A(pix)=ccel;
            I2A(pix0)=0;
            IS6=I2A;
            I22=zeros(size(I2));
            pix1=find(IS6==ccel);
            I2(pix1)=0;
            pix2=unique(I2(:))';
            pix2=pix2(2:end);
            if ccel == 1
                for ity=1:length(pix2)
                  pix4=find(I2==pix2(ity));
                  I22(pix4)=ity;
                end
                TETC{1,im_no}=IS6;
            else
                if ~isempty(pix2)
                    for ity=1:length(pix2)
                      pix4=find(I2==pix2(ity));
                      I22(pix4)=ity;
                    end
                else
                    I22 = I2;
                end
                IS61 = TETC{1,im_no};
                IS61(pix) = ccel;
                TETC{1,im_no}=uint16(IS61);
            end
            TETC{2,im_no}=I22;
        end
    end

    xx = 0;
    for i = rang
        if sum(sum(TETC{2,i})) > 0 && ~isempty(TETC{2,i})
            xx = i;
            break;
        end
    end
    ccel = ccel + 1;
    rang2 = xx:numel(tet_masks);
    disp(xx);

end
ccel = ccel - 1; %number of cells tracked

%% Removing the shock induced points from rang

rang3 = rang;
if ~isempty(shock_period)
    for j = 1:size(shock_period,1)
        for i = shock_period(j,1):shock_period(j,2)
            rang3 = rang3(rang3~=i);
        end
    end
end

%% Removing artifacts - cells that appear once and cells that disapper thresh % of the time or more

all_obj = SR_240222_cal_allob(ccel,TETC,rang);
cell_data = SR_240222_cal_celldata(all_obj,ccel);

k = 1;
cell_artifacts = [];
for iv = 1:ccel
    if cell_data(iv,3) < thresh_perecent*size(rang3,2) || cell_data(iv,5) > thresh
        cell_artifacts(k) = iv;
        k = k +1;
    end
end

all_ccel = 1:ccel;

if ~isempty(cell_artifacts)
    cell_artifacts = unique(cell_artifacts);
    for iv = cell_artifacts
        for its = rang3
            pix = find(TETC{1,its}==iv);
            TETC{1,its}(pix) = 0;
        end
    end
end

%% Retaining and relabeling the new cells

good_cells = setdiff(all_ccel, cell_artifacts);
good_cells = sort(good_cells);

for iv=1:size(good_cells,2)
    for its=rang3
        pix = find(TETC{1,its}==good_cells(iv));
        TETC{1,its}(pix) = iv;
    end
end

% size(good_cells,2) - no. of cells remaining 

%% Correcting the SpoSeg track masks or filling the empty spaces between the first and last appearance
% + Removing artifacts

all_obj1 = SR_240222_cal_allob(size(good_cells,2),TETC,rang);

% for iv = 1:size(all_obj1,1)
%     app = find(all_obj1(iv,:) > 0,2, 'last');
%     x = app(2)-app(1);
%     if x > thresh_remove_last_mask
%         pix = find(TETC{1,app(2)}==iv);
%         TETC{1,app(2)}(pix) = 0;
%         all_obj1(iv,app(2)) = sum(sum(logical(TETC{1,app(2)} == iv)));
%     end
% end

cell_data1 = SR_240222_cal_celldata(all_obj1,size(good_cells,2));


for iv=1:size(good_cells,2)
    for its=cell_data1(iv,1)+1 : cell_data1(iv,2)-1
        if all_obj1(iv,its) == 0
            prev = find(all_obj1(iv,1:its-1) > 0,1, 'last');
            % next = find(all_obj1(iv,its+1:end) > 0,1, 'first')+its;
            % current = its;
            all_obj1(iv,its) = all_obj1(iv,prev);
            pix = find(TETC{1,prev}==iv);
            TETC{1,its}(pix) = iv;
        end
    end
end


%%  cell array that contains the fully tracked TetSeg masks

TETmasks=cell(1,size(TETC,2));
for ita=1:size(TETC,2)
TETmasks{1,ita}=TETC{1,ita};
end

%% Calculate the size of tetrads

TET_obj = size(good_cells,2);
all_obj_final = SR_240222_cal_allob(TET_obj,TETmasks,1:numel(TETmasks));

TET_Size = all_obj_final;

%% Calculate first detection andlast detection of tetrads

TET_exists = zeros(2,TET_obj);
for iv=1:TET_obj
    TET_exists(1,iv) = find(TET_Size(iv,:) > 0,1, 'first'); %1st occurance
    TET_exists(2,iv) = find(TET_Size(iv,:) > 0,1, 'last'); %last occurance
end

%%

tet_masks_exists_tp = rang3;

%%
save([sav_path pos '_TET_Track.mat'],'start','TET_Size','TET_obj','TET_exists','TETmasks','shock_period','thresh','thresh_next_cell','thresh_perecent','tet_masks_exists_tp','-v7.3');   


