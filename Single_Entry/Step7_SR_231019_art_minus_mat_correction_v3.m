
clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

%% Removes the mating events from sposeg tracks based on the overlapping tracked indices 

mat_track_path = [sav_path pos '_MAT_16_18_Track1'];
if ~isempty(dir([mat_track_path, '*']))
    file_list = dir([mat_track_path, '*']);

    mat = load([sav_path file_list(1).name]);

    art_track_path = [sav_path pos '_ART_Track'];
    if ~isempty(dir([art_track_path, '*']))
        file_list = dir([art_track_path, '*']);
        art = load([sav_path file_list(1).name]);

        Mask3 = art.Mask3;
        Matmasks = mat.Matmasks;
        
        for iv = 1:mat.no_obj
            indx_remov = [];
            final_indx_remov = [];
            for its = mat.cell_data(iv,1):mat.cell_data(iv,2) % check for 10 time points
                M = Matmasks{1,its}; %figure;imagesc(M)
                M0 = uint16(M == iv); %figure;imagesc(M0)
                A = Mask3{1,its}; %figure;imagesc(A)
                M1 = imbinarize(M0); %figure;imagesc(M1)
                M2 = bwmorph(M1,'thin',30); %figure(1);imagesc(M2) %imerode(M1,ones(9))
                M3 = A.*M2; %figure(1);imagesc(M3);pause;
                indx = unique(A(M3~=0))';
                if ~isempty(indx)
                    for itt2 = indx
                        if sum(sum(M3==itt2)) > 5
                            indx_remov = [indx_remov itt2];
                        end
                    end
                end
            end
            if ~isempty(indx_remov)
                indx_remov_inter = unique(indx_remov);
                final_indx_remov = unique(indx_remov);
                for itt1 = indx_remov_inter
                    dist_data = -1.*ones(1,numel(Mask3));
                    for its1 = mat.cell_data(iv,1):art.cell_exists(2,itt1)
                        if its1 >= art.cell_exists(1,itt1)
                            M6 = (Mask3{1,its1}==itt1);
                            M7 = (Matmasks{1,its1}==iv);
                            dist_data(1,its1) = sum(sum((M6.*M7)))/sum(sum(M6));
                        end
                    end
                    if ~isempty(find(dist_data(1,:)~=-1,1))
                        first_ov = find(dist_data(1,:)~=-1,1,'first');
                        last_ov = find(dist_data(1,:)~=-1,1,'last');
                        val_avg = median(dist_data(1,first_ov:last_ov));%round(median(dist_data(1,first_ov:last_ov)),1);
                        if val_avg <= 0.4
                            final_indx_remov = setdiff(final_indx_remov,itt1);
                        end
                        % figure(1);plot(dist_data(1,first_ov:last_ov));title(['iv = ' num2str(itt1) ' val = ' num2str(val_avg)]);pause;
                    end
                end
                for its = mat.cell_data(iv,1):numel(Mask3)
                    for itt = final_indx_remov
                        Mask3{1,its}(Mask3{1,its}==itt) = 0;
                    end
                end
            end
            iv
        end
    
        shock_period = mat.shock_period;
        no_obj = art.no_obj;
        ccell2 = art.ccell2;
        cell_exists = art.cell_exists;
        im_no = art.im_no;
    end
    
else
    art_track_path = [sav_path pos '_ART_Track'];
    file_list = dir([art_track_path, '*']);
    filename_old = file_list(1).name;
    filename_new = replace(filename_old,'_ART_Track','_ART_Track1');
    copyfile([sav_path filename_old],[sav_path filename_new]);
end
        

save([sav_path pos '_ART_Track1.mat'],"no_obj","shock_period","Mask3","im_no","ccell2","cell_exists", "-v7.3")


