%% CORRECTING ProSeg Masks USING SpoSeg Tracks

clc;
clear;
close all;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\']; % path to the segmented tet masks
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';
shock_period = [122,134];

path_dir = dir(fullfile(path, '*_ART_masks.tif'));
Art_MT = cell(1,numel(path_dir));
for its = 1:numel(path_dir)
    Art_MT{1,its} = uint16(imread(fullfile(path,path_dir(its).name)));
end
    
tet_track_path = [sav_path pos '_TET_Track.mat'];
if ~isempty(dir([tet_track_path, '*']))
    file_list = dir([tet_track_path, '*']);
    tet = load([sav_path file_list(1).name]); % Load the tracked SpoSeg masks
    shock_period = tet.shock_period;

    for iv = 1:tet.TET_obj % Loop through the total number of TETs detected
        if tet.TET_exists(2,iv) >= shock_period(1,1)-1 % Conditional loop since we want to focus on only the time points before the shock
            tp_end = shock_period(1,2);
        else
            tp_end = tet.TET_exists(2,iv);
        end
        for its = tet.TET_exists(1,iv):tp_end
            A1 = double(Art_MT{its}); % figure;imagesc(A1)
            if its >= shock_period(1,1) && its <= shock_period(1,2)
                T1 = double(tet.TETmasks{1,shock_period(1,1)-1}==iv); % figure;imagesc(T1)
                thresh = 0.6;
            else
                T1 = double(tet.TETmasks{1,its}==iv);
                thresh = 0.95;
            end
            T1 = double(imresize(T1,size(A1),'nearest'));
            Im1=imbinarize(T1); % figure(3);imagesc(Im1)
            Im2=imerode(Im1,ones(9));% figure;imagesc(Im2)
            Im3=A1.*Im2;  %figure;imagesc(Im3)
    
            pix11 = [];
            pix1=unique(A1(Im3~=0))';
            for it2 = pix1
                r1 = sum(sum(Im3==it2))/sum(sum(logical(Im3)));
                if r1 > 0.2
                    pix11 = [pix11 it2];
                end
            end
    
            if size(pix11,2) == 1
            r = sum(sum(A1==pix11))/sum(sum(T1));
            if r > thresh
            else
                Art_MT{1,its}(A1(:)==pix11) = 0;
                Art_MT{1,its}(T1(:)==1) = max(Art_MT{1,its}(:))+1;
            end                
            elseif isempty(pix11) %!!!!!!!!!!!!! ART CANNOT BE MODIFIED FOR MULTIPLE ENTRY PAY ATTENTION 
                    Art_MT{1,its}(T1(:)==1) = max(Art_MT{1,its}(:))+1;
            else
                for it2=pix11
                    Art_MT{1,its}(A1==it2)=0; % Removing all the indexes from the ART mask
                end
                Art_MT{1,its}(T1==1)=max(Art_MT{1,its}(:))+1; %Replacing the ART mask with TET mask
            end
        end
    end
    
    
    for iv = 1:tet.TET_obj % Loop through the total number of TETs detected
        if tet.TET_exists(2,iv) > shock_period(1,2) && tet.TET_exists(1,iv) < shock_period(1,1) % Conditional loop since we want to focus on only the time points before the shock
            % figure(1);plot(tet.TET_Size(iv,shock_period(1,2)+1:tet.TET_exists(2,iv)));pause;
            s1 = sum(sum(tet.TETmasks{1,shock_period(1,2)+1}==iv));
            for its = shock_period(1,2)+1:tet.TET_exists(2,iv)
                A1 = double(Art_MT{its}); %figure(1);imagesc(A1);%pause;
                T1 = double(tet.TETmasks{1,its}==iv); %figure(2);imagesc(T1)
                  
                s2 = sum(sum(tet.TETmasks{1,its}==iv));
                if its == tet.TET_exists(2,iv)
                    s3 = sum(sum(tet.TETmasks{1,its}==iv));
                else
                    s3 = sum(sum(tet.TETmasks{1,its+1}==iv)); %future mask
                end
                if s2 < s1 - 0.1*s1
                    if s3 > s2 + 0.1*s2
                        T1 = double(tet.TETmasks{1,its-1}==iv);
                    else
                    % disp(iv);
                    % disp(its);
                    break;
                    end
                end
                s1 = s2;
                T1 = double(imresize(T1,size(A1),'nearest'));
                Im1=imbinarize(T1); % figure(3);imagesc(Im1)
                Im2=imerode(Im1,ones(9));% figure;imagesc(Im2)
                Im3=A1.*Im2; %figure;imagesc(Im3)
                
                pix11 = [];
                pix1=unique(A1(Im3~=0))';
                for it2 = pix1
                    r1 = sum(sum(Im3==it2))/sum(sum(logical(Im3)));
                    if r1 > 0.2
                        pix11 = [pix11 it2];
                    end
                end
    
                if size(pix11,2) == 1
                r = sum(sum(A1==pix11))/sum(sum(T1));
                if r > thresh
                else
                    Art_MT{1,its}(A1(:)==pix11) = 0;
                    Art_MT{1,its}(T1(:)==1) = max(Art_MT{1,its}(:))+1;
                end                
                elseif isempty(pix11)
                        Art_MT{1,its}(T1(:)==1) = max(Art_MT{1,its}(:))+1;
                else
                    for it2=pix11
                        Art_MT{1,its}(A1==it2)=0; % Removing all the indexes from the ART mask
                    end
                    Art_MT{1,its}(T1==1)=max(Art_MT{1,its}(:))+1; %Replacing the ART mask with TET mask
                end
                % figure(1);imagesc(Art_MT{1,its});title(['Masks - tp = ' num2str(its)]); pause;
            end
        end
    end

    save([sav_path pos '_ART_Masks.mat'] ,"Art_MT","shock_period","-v7.3");
else
    % shock_period = [];
    save([sav_path pos '_ART_Masks.mat'] ,"Art_MT","shock_period","-v7.3");
end
        
