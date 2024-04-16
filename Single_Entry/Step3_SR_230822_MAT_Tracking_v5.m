clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\']; % path to the segmented mat masks
sav_path = 'E:\SR_Tracking\toy_data\Tracks\'; % Path to save Tracks
shock_period = [122,134];

path_dir=dir(fullfile(path,'img_*_Ph3_000_MAT16_18_masks.tif'));
fileNames = {path_dir.name};
fileNumbers = zeros(size(fileNames));

for i = 1:numel(fileNames)
    fileNumbers(i) = sscanf(fileNames{i}, 'img_%d_Ph3_000_MAT16_18_masks.tif');
end

[sortedNumbers, sortedIndices] = sort(fileNumbers);

mat_masks_path = {};
for i = sortedIndices
    img_path = fullfile(path, path_dir(i).name);
    mat_masks_path = [mat_masks_path, img_path];
end

mat_masks = {};
for i = 1:numel(fileNames)
    % +1 because the images start from 0
    mat_masks{sortedNumbers(i)+1} = imread(mat_masks_path{i});
end

% Save all the valid mat masks in a variable called mat_masks
for i = min(sortedNumbers)+1:size(mat_masks,2)
    if isempty(mat_masks{i})
        mat_masks{i} = uint16(zeros(size(mat_masks{min(sortedNumbers)+1})));
    end
end

%% Remove shock induced timepoints

mat_masks_original = mat_masks;

if ~isempty(shock_period)
    for j = 1:size(shock_period,1)
        %for i = shock_period(j,1):shock_period(j,2)
        for i = 1:shock_period(j,2)
            mat_masks{i} = [];    
        end
    end
end

start = 0;
for its=1:size(mat_masks,2)
    M = uint16(mat_masks{its});
    if sum(sum(M)) > 0
        start = its;
        break;
    end
end

    %% Tracking all the detections

if start ~=0
    rang = start:numel(mat_masks);
    I2 = mat_masks{start};
    A = zeros(size(mat_masks{start}));
else
    rang = 1:numel(mat_masks);
    I2 = mat_masks{1};
    A = zeros(size(mat_masks{1}));
end

IS6=zeros(size(I2));
MATC = cell(2,numel(mat_masks));
xx = start;
rang2 = rang;
ccel = 1;

while xx ~= 0
    for im_no=rang2
        if ccel == 1
            I2 = mat_masks{im_no}; % read the detected mat masks
        else
            I2 = MATC{2,im_no};
        end

        if isempty(I2) % figure;imagesc(I2)
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
            % overlap from previous time point
            I3B =uint16(I3A).*uint16(I2A);% figure;imagesc(uint16(I3A));figure;imagesc(uint16(I2A))
            ind = mode(I3B(I3B~=0));
            if ind==0 && ccel == 1
                MATC{1,im_no}=I3B;
                MATC{2,im_no}=I2A;
                continue;
            elseif ind==0 && ccel~=1
                continue;
            end
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
                MATC{1,im_no}=IS6;
            else
                if ~isempty(pix2)
                    for ity=1:length(pix2)
                      pix4=find(I2==pix2(ity));
                      I22(pix4)=ity;
                    end
                else
                    I22 = I2;
                end
                IS61 = MATC{1,im_no};
                IS61(pix) = ccel;
                MATC{1,im_no}=uint16(IS61);
            end
            MATC{2,im_no}=I22;
        end
    end
    xx = 0;
    for i = rang
        if sum(sum(MATC{2,i})) > 0 && ~isempty(MATC{2,i})
            xx = i;
            break;
        end
    end
    % disp(ccel);
    disp(xx);
    ccel = ccel + 1;
    rang2 = xx:numel(mat_masks);
end
ccel = ccel - 1;

%% Removing the shock induced points from rang

rang3 = rang;
if ~isempty(shock_period)
    for j = 1:size(shock_period,1)
       % for i = shock_period(j,1):shock_period(j,2) %%!!!!!!!!!!!!!!!!!
        for i = 1:shock_period(j,2)    
            rang3 = rang3(rang3~=i);
        end
    end
end

%% Correction Code

all_obj = SR_240222_cal_allob(ccel,MATC(1,:),rang);
cell_data = SR_240222_cal_celldata(all_obj,ccel);


for iv = 1:ccel
    if ~isempty(find(all_obj(iv,min(rang):shock_period(size(shock_period,1),2))>0,1))
        if all_obj(iv,shock_period(size(shock_period,1),2)+1) ~= 0
            for its = shock_period(size(shock_period,1),2)+1:rang(end)
                if all_obj(iv,its)~= -1
                    pix = find(MATC{1,its}==iv);
                    MATC{1,its}(pix) = 0;
                    all_obj(iv,its) = sum(sum(logical(MATC{1,its} == iv)));
                end
            end
        end
    end
end

cell_data = SR_240222_cal_celldata(all_obj,ccel);

k = 1;
cell_artifacts = [];
for iv = 1:ccel
    if cell_data(iv,3) == 1 || cell_data(iv,5) > 80
        cell_artifacts(k) = iv;
        k = k +1;
    end
end

all_ccel = 1:ccel;

if ~isempty(cell_artifacts)
    cell_artifacts = unique(cell_artifacts);
    for iv = cell_artifacts
        for its = rang3
            pix = find(MATC{1,its}==iv);
            MATC{1,its}(pix) = 0;
        end
    end
end

good_cells = setdiff(all_ccel, cell_artifacts);
good_cells = sort(good_cells);

for iv=1:size(good_cells,2)
    for its=rang3
        pix = find(MATC{1,its}==good_cells(iv));
        MATC{1,its}(pix) = iv;
    end
end

ccel = size(good_cells,2);
all_obj = SR_240222_cal_allob(ccel,MATC(1,:),rang);
cell_data = SR_240222_cal_celldata(all_obj,ccel);

for iv=1:ccel
    tp_data{iv,1} = diff(find(all_obj(iv,:)>0));
    tp_data{iv,2} = find(all_obj(iv,:)>0);
    a = find(tp_data{iv}>10);
    if ~isempty(a)
        if a == size(tp_data{iv},2)
            pix = find(MATC{1,tp_data{iv,2}(a+1)}==iv);
            MATC{1,its}(pix) = 0;
        else
            for its = find(all_obj(iv,:) > 0,1, 'first'):tp_data{iv,2}(a+1)-1
                pix = find(MATC{1,its}==iv);
                MATC{1,its}(pix) = 0;
            end
        end
    end
end

for iv=1:ccel
    for its=find(all_obj(iv,:)>0,1,'first')+1 : find(all_obj(iv,:)>0,1,'last')-1
        if all_obj(iv,its) == 0
            prev = find(all_obj(iv,1:its-1) > 0,1, 'last');
            % next = find(all_obj(iv,its+1:end) > 0,1, 'first')+its;
            % current = its;
            all_obj(iv,its) = all_obj(iv,prev);
            pix = find(MATC{1,prev}==iv);
            MATC{1,its}(pix) = iv;
        end
    end
end

all_obj = SR_240222_cal_allob(ccel,MATC(1,:),rang);
cell_data = SR_240222_cal_celldata(all_obj,ccel);

no_obj = ccel;
Matmasks=cell(1,size(rang,2));
for ita=rang
    Matmasks{1,ita}=MATC{1,ita};
end

save([sav_path pos '_MAT_16_18_Track.mat'],"Matmasks","no_obj","all_obj","cell_data","rang","rang3","shock_period","mat_masks_original","start",'-v7.3');



