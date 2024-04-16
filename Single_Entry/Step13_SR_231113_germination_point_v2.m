clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\'];
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';

tet_id_path = [sav_path, pos, '_TET_ID_art_track.mat'];
file_list = dir([tet_id_path, '*']);
load([sav_path file_list(1).name]); % load tet IDs

art_track1_path = [sav_path, pos, '_ART_Track1_DS.mat'];
file_list = dir([art_track1_path, '*']);
load([sav_path file_list(1).name]); % load art track

tet_track_path = [sav_path, pos, '_TET_Track_DS.mat'];
file_list = dir([tet_track_path, '*']);
tet = load([sav_path file_list(1).name]); % load tet track

int = tet.shock_period(1,2)+1:tet.shock_period(1,2)+35; % choose time frame during rich medium where cells will likely germinate

desc_path = [sav_path, pos, '_final_descendants.mat'];
file_list = dir([desc_path, '*']);
desc = load([sav_path file_list(1).name]); % load descendants information
alive_tets = desc.alive_tets;

germination_point = zeros(1,size(TET_ID,2));
for iv = alive_tets%1:size(TET_ID,2) % loop through the tet cells that are alive
    if TET_ID(1,iv) ~= -1 % if it is a tet that succesfully germinates and has descendants
        y = all_ob(TET_ID(1,iv),:);
        A = all_ob(TET_ID(1,iv),int);
        [v,i] = max(A);
        k = i+tet.shock_period(1,2);
        B = y(k:k+5);
        C = diff(B); % the drastic drop in cell size after germination defines the germination point
        [iP, iL] = max(-C);
        k1 = k + iL - 1;
        germination_point(1,iv) = k1; %is value = 0, then tet is dead
    end
end

save([sav_path pos '_germination_point'],"germination_point","alive_tets");
