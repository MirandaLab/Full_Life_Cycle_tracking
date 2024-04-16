%% THIS IS FOR ART TRACK 1 TRACKS

clc;
clear;
% close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

% Load TET Tracks
tet_track_path = [sav_path pos '_TET_Track_DS'];
file_list = dir([tet_track_path, '*']);
load([sav_path file_list(1).name]);  % load tet track

tet_track_path = [sav_path pos '_TET_ID'];
file_list = dir([tet_track_path, '*']);
load([sav_path file_list(1).name]); % load tet IDs

% Load ART Tracks
art_track_path = [sav_path pos '_ART_Track_DS'];
file_list = dir([art_track_path, '*']);
load([sav_path file_list(1).name]); %  load art track

masks = TETmasks;
ART = Mask3;
shock_period = shock_period;
for i = 1:numel(masks)
    if ~isempty(masks{i})
    masks{i} = imresize(masks{i},size(ART{i}),"nearest");
    end
end

start = shock_period(1,2)+1;
tp_end = numel(ART); %190;
% nutr_stress_tp = numel(ART);%190;
% intx2 = start:nutr_stress_tp;
int = start:tp_end;

%% To determine if the tet cells are dead or alive after shock period

size_var_tet = zeros(1,TET_obj);
dead_tets = zeros(1,TET_obj);
for iv = 1:TET_obj
     if TET_ID(1,iv) ~= -1
        A = all_ob(TET_ID(1,iv),shock_period(1,2)+1:end)'; %size is calculated only from shock end +1 to the end
        [LT,ST,R] = trenddecomp(A); % plot([A LT ST R]);
        size_var_tet(1,iv) = var(ST(:,1)); %Calculate the variance of the seasonal trend
        if var(ST(:,1)) < 1000 %thresh
            dead_tets(iv) = 1;
        else
            dead_tets(iv) = 0;
        end
     else
         size_var_tet(1,iv) = -10.^5; %if they have a TET_ID of -1
         dead_tets(iv) = -10.^5; %if they have a TET_ID of -1
     end
end

%% Finding the first time when tets germinate and new cells begin to show up

begin = 0;
for its = int
    % its = min(int)
    if its == max(int)
        begin = 0;
        break;
    end
    A1 = ART{1,its}; %imagesc(A1)
    A2 = ART{1,its+1}; %imagesc(A2)
    A3 = uint16(logical(A1)).*uint16(A2); %imagesc(A3)
    indx_ori = unique(A1(A1~=0)); % previous mask
    indx_new = unique(A2(A2~=0)); % present mask
    vals = setdiff(indx_new,indx_ori);
    if ~isempty(vals)
        begin = its+1; % where new cells begin to emerge
        break;
    end
end

%% Dividing the FOV into regions that belong to certain tetrads using the watershed algorithm

if begin ~= 0

    int1 = begin:numel(ART)-1;
    new_indx = cell(1,numel(ART));
    
    for its = int1
        A1 = ART{1,its-1}; %imagesc(A1)
        A2 = ART{1,its}; %imagesc(A2)
        indx_ori = unique(A1(A1~=0)); % previous mask
        indx_new = unique(A2(A2~=0)); % present mask
        vals = setdiff(indx_new,indx_ori);
        new_indx{its} = vals';
    end
    kka = 0; %Add this after TK_OAM_SKI_whi5_230330 pos = 59 Pos8_1
    new_indx_new = {};
    for ii1 = 1:numel(new_indx)
        if ~isempty(new_indx{1,ii1})
            kka = kka+1;
            new_indx_new{kka} = new_indx{1,ii1};
        end
    end
    
    new_born = unique(cell2mat(new_indx_new));
    
    I2 = zeros(size(ART{start})); % at this point both art and tet masks have the same size
    for ccell = 1:TET_obj
        if TET_ID(1,ccell) ~= -1
            if dead_tets(ccell) == 0
                if TET_exists(ccell,2) >= shock_period(1,2)+1
                    stats = regionprops(masks{shock_period(1,2)+1}==ccell,'Centroid');
                else
                    stats = regionprops(masks{TET_exists(ccell,2)}==ccell,'Centroid');
                end
        cent = round(stats(1).Centroid); %[x,y]
        I2(cent(2),cent(1)) = 1;
        % I2(cent(2)-5:cent(2)+5,cent(1)-5:cent(1)+5) = 0;
            end
        end
    end
    
    I21 = bwmorph(I2,'thicken',ones(9,9));
    I4 = bwdist(I21); %figure;imagesc(I4)
    I3 = watershed(I4,4); %figure;imagesc(I3) % Creating the boundaries
    
%% Checking which cell from art tracks belongs to which region created using watershed

    region = [];
    amt = []; % 
    k = 0;
    for iv = 1:no_obj
        I12 = zeros(size(uint16(ART{1,start})));
        kx = 0;
        for its = int1
            I11 = logical(ART{1,its}==iv);
            if sum(I11(:)) > 0
                kx = kx+1;
                if kx >= 1 && kx <= 2 %kx==1
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
            region(k,1) = iv; % cell number in art tracks
            region(k,2) = pix(ind); %region
        end
    end
    
    unique_regions = unique(region(:, 2)); %% All cells in artilife that belong to a particular sector
    cell_arrays = cell(length(unique_regions), 1);
    for i = 1:length(unique_regions)
        current_region = unique_regions(i);
        cells_in_current_region = region(region(:, 2) == current_region,1);
        cell_arrays{i} = uint16(cells_in_current_region); % 
    end
    
%% Saving the tet id, tet regions and possible descendants

    cell_re = {}; % !!!!!!! GO THROUGH
    TET_ind = zeros(numel(unique_regions),3);
    common_indices = cell(1,TET_obj);
    amt_1 = [];
    for iv = 1:TET_obj % loop through tet indicies
        if TET_ID(1,iv) ~= -1
            if dead_tets(iv) == 0
        
        T1 = masks{1,TET_exists(iv,1)}==iv; % figure;imagesc(T1);
        T2 = uint16(I3).*uint16(T1); % figure;imagesc(T2);
        % % T3 = uint16(ART{1,TET_exists(iv,1)}).*uint16(T1); %figure;imagesc(T3);
        % T3 = uint16(ART{1,shock_period(1,2)+1}).*uint16(T1);
        % if sum(T3(:)) == 0
        %     T3 = uint16(ART{1,shock_period(1,1)-1}).*uint16(T1);
        %     if sum(T3(:)) == 0
        %         T3 = uint16(ART{1,TET_exists(iv,2)}).*uint16(T1);
        %     end
        % end
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
    
        common_indices{iv} = intersect(new_born,cell_arrays{TET_ind(iv,2)});
        common_indices{iv} = [common_indices{iv}; TET_ind(iv,3)]; % if we want to include the tet index as well
            
            else
                TET_ind(iv,1) = iv; % TET number
                TET_ind(iv,3) = TET_ID(1,iv);
                common_indices{iv} = [common_indices{iv}; TET_ind(iv,3)];
            end
        end
    end
%% This part of the code works on identifying incorrectly associated descedndats, the code removes the association when incorrect and assigns it to the right tet based on connectivity

    alive_tets = [];
    for iv = 1:TET_obj
        if dead_tets(iv) == 0
            alive_tets = [alive_tets iv];
        end
    end
    
    k = 0;
    kk = 0;
    need_remov = [];
    works = [];
    for iv = alive_tets
        common_indices1 = common_indices{iv}';
        for ittx1 = common_indices{iv}' %intx2!!!!!!!!!!!!!!
            if ~isempty(find(TET_ID(1,:)==ittx1,1)) %ignore TET trads
                continue
            end
            its = cell_exists(1,ittx1);
            M = uint16(ART{1,its}); %imagesc(M) %!!!!(1)
            I_s0= uint16(zeros(size(M))); %!!! check doubles(2)
            I_s2 = uint16(zeros(size(M)));
            for itx1 = 1:size(need_remov,2)
                common_indices1 = setdiff(common_indices1,need_remov(1,itx1));
            end
            for it=common_indices1
                I_s2=(M==it); % imagesc(I_s1)
                I_s2=bwmorph(I_s2,'thicken',5);%imagesc(I_s0); % try bwmorph with gpuArray
                I_s0=I_s0+uint16(I_s2);
            end
            IA1 = imbinarize(I_s0);
            IA2 = bwlabel(IA1); %figure(1);imagesc(IA2);title(['iv = ' num2str(ittx1)]);pause; % gpu
            max(IA2(:))
            if max(IA2(:)) > 1
                out = unique(IA2(:))';
                sizes_occup = []; 
                for itt1 = out
                    if (itt1~=0)
                    a = sum(sum(IA2 ==itt1));
                    sizes_occup = [sizes_occup a]; 
                    end
                end
                xx1 = 0;
                for itt1 = out
                    if (itt1~=0)
                        IA3 = uint16(imerode(IA2==itt1,strel('disk',5)));
                        IA4 = IA3.*uint16(ART{1,its});
                        pixx5 = unique(IA4);
                        if numel(pixx5)>1
                            if numel(pixx5)>2 || ~isempty(find(TET_ID == pixx5(2), 1))
                                if sum(sum(IA2==itt1)) == max(sizes_occup)
                                    xx1 = itt1;
                                end
                            end
                        end
                    end
                end
                AAB = 0;
                for itt2 = out
                    if itt2 ~= 0
                        AB1 = uint16(IA2);
                        AB2 = uint16(bwmorph(M == ittx1,'thin',5));
                        AB3 = AB1.*AB2;
                        pixab = unique(AB3)';
                        if numel(pixab) > 1
                            if pixab(2) == xx1
                                continue;
                            else
                                AAB = 1;
                            end
                        elseif numel(pixab) == 1
                            continue;
                        else
                                disp(ittx1);
                                disp(iv);
                        end
                    end
                end
                if AAB == 1
                    for itx = alive_tets
                        if itx ~= iv
                            M = uint16(ART{1,its});
                            I_s0= uint16(zeros(size(M)));
                            I_s2 = uint16(zeros(size(M)));
                            for it=[common_indices{itx}',ittx1]
                                I_s2=(M==it); % imagesc(I_s1)
                                I_s2=bwmorph(I_s2,'thicken',3);% imagesc(I_s0)
                                I_s0=I_s0+uint16(I_s2);
                            end
                            IA11 = imbinarize(I_s0);
                            IA21 = bwlabel(IA11); %imagesc(IA21);pause;
                            if max(IA21(:)) == 1
                                k = k + 1;
                                works(1,k) = ittx1;
                                works(2,k) = iv; % old association
                                works(3,k) = itx; % new association
                            else
                                kk = kk+1;
                                need_remov(1,kk) = ittx1;
                                need_remov(2,kk) = iv; % old association
                            end
                        end
                    end
                end
            end
        end
    end    
    
    descendants = common_indices;
    for itx1 = 1:size(works,2)
        descendants{works(2,itx1)} = setdiff(descendants{works(2,itx1)},works(1,itx1));
        descendants{works(3,itx1)} = [descendants{works(3,itx1)};works(1,itx1)];
    end
    for itx1 = 1:size(need_remov,2)
        descendants{need_remov(2,itx1)} = setdiff(descendants{need_remov(2,itx1)},need_remov(1,itx1));
    end
    
    descendants_data = cell(TET_obj,3);
    for iv = 1:TET_obj % loop through tet indicies
        if TET_ID(1,iv) == -1
            descendants_data{iv,1} = iv;
            descendants_data{iv,2} = TET_ID(1,iv);
            descendants_data{iv,3} = -1;
        else
            descendants_data{iv,1} = iv;
            descendants_data{iv,2} = TET_ID(1,iv);
            descendants_data{iv,3} = setdiff(descendants{iv},TET_ID(1,iv));
        end
    end
    
    save([sav_path pos '_descendants_new_art.mat'],"I3","descendants_data","descendants","alive_tets","common_indices","cell_arrays","TET_obj","-v7.3");

end