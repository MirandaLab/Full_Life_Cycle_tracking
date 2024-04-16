% Author: Shreya Ramakanth

clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

art_track_path = [sav_path pos '_ART_Track1_DS'];
file_list = dir([art_track_path, '*']);
art = load([sav_path file_list(1).name]);

mat_track_path = [sav_path pos '_MAT_16_18_Track1_DS'];
file_list = dir([mat_track_path, '*']);
mat = load([sav_path file_list(1).name]);

%% Obtain the gametees indexes that give rise to the mating cell (zygote)

Mask3 = art.Mask3; %Reading Art Tracks from art
MTrack = mat.Matmasks; %Reading Mat Tracks from mat
gamete=zeros(2,mat.no_obj);

for iv = 1:mat.no_obj
    tp_mat_start = mat.cell_data(iv,1); % First appearance of mating event "iv"
    M1 = MTrack{1,tp_mat_start}; % Mat Track at tp_mat_start
    for its = tp_mat_start-1:-1:mat.shock_period(1,2)+1 % Loop through time from 1 tp before mating to one time point after shock
        A1 = double(Mask3{1,its});% figure();imagesc(A1)
        M2 = double(M1 == iv);

        Im1=imbinarize(M2); % figure(3);imagesc(Im1)
        Im2=bwmorph(Im1,'thin',10);% figure;imagesc(Im2)
        Im3=A1.*Im2; % figure;imagesc(Im3)
        pix2=unique(A1(Im3~=0))';
        % figure(1);imagesc(M2);
        % figure(2);imagesc(A1);
        % figure(3);imagesc(Im3);pause;

        if size(pix2,2) == 2 % captures mature mating
            r = sum(sum(Im3==pix2(1,1)))/sum(sum(Im3==pix2(1,2)));
            if ((2/8) <= r) && (r <= (8/2)) %4/6 to 9/6
                gamete(1,iv) = pix2(1,1);
                gamete(2,iv) = pix2(1,2);
                gamete(3,iv) = its;
            end
        end
    end
end

for iv = 1:mat.no_obj
    if gamete(1,iv) == 0 && gamete(2,iv) == 0
        tp_mat_start = mat.cell_data(iv,1); % First appearance of mating event "iv"
        M1 = MTrack{1,tp_mat_start}; % Mat Track at tp_mat_start
        for its = tp_mat_start-1:-1:mat.shock_period(1,2)+1 % Loop through time from 1 tp before mating to one time point after shock
            A1 = double(Mask3{1,its});% figure();imagesc(A1)
            M2 = double(M1 == iv);

            Im1=imbinarize(M2); % figure(3);imagesc(Im1)
            Im2=bwmorph(Im1,'thin',10);% figure;imagesc(Im2)
            Im3=A1.*Im2; % figure;imagesc(Im3)
            pix2=unique(A1(Im3~=0))';
            % figure(1);imagesc(M2);
            % figure(2);imagesc(A1);
            % figure(3);imagesc(Im3);pause;

            if size(pix2,2) == 1 % captures ascus mating
                gamete(1,iv) = pix2(1,1);
                gamete(3,iv) = its;
            end
        end
    end
end
save([sav_path pos '_gametes'],"gamete",'-v7.3');

