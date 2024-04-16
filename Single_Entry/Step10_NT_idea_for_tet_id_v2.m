
clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

art_path = [sav_path pos '_ART_Track_DS'];
art_file_list = dir([art_path, '*']);
load([sav_path art_file_list(1).name]);  % load art track

tet_path = [sav_path pos '_TET_Track_DS'];
tet_file_list = dir([tet_path, '*']);
load([sav_path tet_file_list(1).name]);  % load tet track


%% Finding the corresponding cell number of tet cells ( in tet tracks) in art atracks

TET_ID = zeros(1,TET_obj);
for iv = 1:TET_obj
    its = TET_exists(iv,1);

    if its>shock_period(1,2)
        TET_ID(1,iv) = -1;
    else 

    A = double(Mask3{its}); % figure;imagesc(A)
    T = double(imresize(TETmasks{its},size(A),'nearest')); % figure;imagesc(T)
    
    T1 = double(T == iv); % figure;imagesc(T1)
    
    Im1=imbinarize(T1); % figure(3);imagesc(Im1)
    Im2=imerode(Im1,ones(9));% figure;imagesc(Im2)
    Im3=A.*Im2; % figure;imagesc(Im3)
    
    % find indexes to erase
    pix1=unique(A(Im3~=0));
    TET_ID(1,iv) = pix1;

    end 
end

name1=[sav_path pos '_TET_ID_art_track'];
save(name1, 'TET_ID');


