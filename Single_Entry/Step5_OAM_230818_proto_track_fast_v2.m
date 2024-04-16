%% Track ProSeg Cells

clc;
clear;
close all;

Arti_v=11;
cell_prob=0.5;%0.5;
flow_threshold=0.9;%0.9;

pos = 'Pos0_2';
path = ['E:\SR_Tracking\toy_data\' pos '\']; % path to the segmented tet masks
sav_path = 'E:\SR_Tracking\toy_data\Tracks\';
art_mask_path = [sav_path pos '_ART_Masks'];
file_list = dir([art_mask_path, '*']);
load([sav_path file_list(1).name]);
Masks3 = Art_MT;

im_no1 = 1;
im_no=size(Masks3,2);% !!!!!!!!!!%size(Masks3,1)-1; %last time point (starting from 0)/total number of timepoints tracked
mm=[1 im_no]; % time points to track

%% load the first mask that begins the indexing for all the cells; IS1 is updated to most recently processed tracked mask at the end of it0
IS1 = uint16(Masks3{1,im_no1});%start tracking at first timepoint % figure; imagesc(IS1);
IS1 = OAM_230919_remove_artif(IS1);% remove artifacts and start tracking at first timepoint % figure; imagesc(IS1);
masks=zeros(size(IS1,1),size(IS1,2),im_no); % contains the re-labeled masks accoridng to labels in the last tp mask
masks(:,:,im_no1) = IS1; %first timepoint defines indexing; IS1 is first segmentation output

%% Allocate a mask for cells where there's a gap in the segmentation; IblankG is updated within the loops it1 & itG
IblankG=zeros(size(IS1)); % title('IblankG'); colorbar;
% tic
for it0=mm(1):mm(2)%im_no1:im_no %206206 % % notice IS1 will be updated in the loop
    
disp(['it0=' num2str(it0)])
% it0=98;
%% load the future cellpose mask, IS2: IS2 is the current image being re-indexed for tracking
IS2=uint16(Masks3{1,it0});% figure;imagesc(IS2==1); title('IS2');figure;imagesc(IS2);
IS2 = OAM_230919_remove_artif(IS2);% remove artifacts

IS2C = IS2;% figure;imagesc(IS2C); title('IS2C') <--- a copy of IS2, gets updated in it1
[IS1B]=OAM_231216_bina(IS1);

IS3 = uint16(IS1B).*IS2;% figure;imagesc(IS3) %past superimposed onto future; updated in it1
% tic
tr_cells=unique(IS1(IS1~=0))'; % tr_cells=tr_cells(2:end)';%the tracked cells present in the past mask, IS1.
% toc
gap_cells=unique(IblankG(IblankG~=0))'; %gap_cells=gap_cells(2:end)'; %the tracked cells that had a gap in their segmentation; were not detected in IS1. %figure; imagesc(IblankG);
cells_tr=[tr_cells gap_cells]; %all the cells that have been tracked up to this tp for the position

%% Allocate space for the re-indexed IS2 according to tracking
Iblank0=zeros(size(IS1));% figure;imagesc(Iblank0);colorbar;

%% Go to the previously tracked cells and find corresponding index in current tp being processed, IS2 -> Iblank0: mask of previously tracked cells with new position in IS2
if sum(cells_tr(:))~=0 % this is required in case the mask goes blanck because cells mate inmediately during germination
for it1=sort(cells_tr) % <---cells are proccessed in order according to birth/appearance 
% it1 =24;
IS5=(IS1==it1);% figure;imagesc(IS5) % go to past mask,IS1, to look for the cell
IS6A=uint16(bwmorph(IS5,'thin')).*IS3;% figure;imagesc(IS6A);

if sum(IS5(:))==0 %if the cell was missing in past mask; look at the gaps in segmentation otherwise continue to look at the past mask
IS5=(IblankG==it1);
IS6A=uint16(bwmorph(IS5,'thin')).*IS2C; 
IblankG(IblankG==it1)=0; %remove the cell from segmentation gap mask - it'll be in the updated past mask for next round of processing
else
end

% Find the tracked cell's corresponding index in IS2, update IS3 and IS2C to avoid overwriting cells 
if sum(IS6A(:))~=0
IS2ind=mode(IS6A(IS6A~=0));
Iblank0(IS2==IS2ind)=it1;% figure; imagesc(Iblank0)
IS3(IS3==IS2ind)=0;% figure; imagesc(IS3)
IS2C(IS2==IS2ind)=0;% figure; imagesc(IS2C)
else
end
end
% %figure; imagesc(Iblank0);
% figure; imagesc(IS2);
% figure; imagesc(IS3);

%% define cells with segmentation gap, update IblankG, the segmentation gap mask
seg_gap=setdiff(tr_cells,(unique(Iblank0)')); %cells in the past mask (IS1), that were not found in IS2 

%IblankG=zeros(size(IS1));
if ~isempty(seg_gap)
for itG=seg_gap
    IblankG(IS1==itG)=itG;
end
else
end
%figure; imagesc(IblankG);

%% define cells that were not relabelled in IS2; these are the buds and new cells entering the frame
Iblank0B=Iblank0; Iblank0B(Iblank0~=0)=1;% figure; imagesc(Iblank0B);
ISB=IS2.*uint16(~(Iblank0B)); %figure; imagesc(~Iblank0B);
%ISBE=IS1.*uint16(~(Iblank0B));%figure; imagesc(IS2);
%   figure; imagesc(ISB);

%% Add new cells to the mask with a new index Iblank0->Iblank, Iblank0 with new cells added
newcells=unique(ISB(ISB~=0))'; %newcells=newcells(2:end)';
Iblank=Iblank0;
A=1;
% figure; imagesc(IS2);
if ~isempty(newcells)
for it2=newcells %figure;imagesc(IS6B)
    Iblank(IS2==it2)=max(cells_tr)+A; % figure;imagesc(Iblank); create new index that hasn't been present in tracking
    A=A+1;
end
else
end

masks(:,:,it0)=uint16(Iblank); %<---convert tracked mask to uint16 and store
IS1=masks(:,:,it0);% IS1, past mask, is updated for next iteration of it0

else

masks(:,:,it0)=IS2; %
IS1=IS2;% IS1, past mask, 

end

end
% toc


%% tracks as as tensor
% tic
im_no=size((masks),3); % last time point to segment   size(masks2)
ccell2=unique(masks(masks~=0))';% !!!!!!!!!!!!! %unique(masks{1,im_no})'; % figure;imagesc(masks{1,im_no})

Mask2=zeros(size(masks)); % 4246 in 3h
    % re-assign ID
    % pix2=cell(1,size(ccell2,2));
    mak=gpuArray(masks);
    for itt3 = 1:size(ccell2,2) % cells
        pix3=find(mak==ccell2(itt3)); % find(cell2mat(pix2)~=0)
        Mask2(pix3)=itt3; 
        itt3
    end


 % toc

%% get cell presence
Mask3=Mask2;
numbM=im_no;
obj=unique(Mask3); %figure; imagesc((cell2mat(masks))); <--- look at timeseries of masks (view with caution) 
no_obj1=max(obj); 
A=1;

 tic % Elapsed time is 9.671897 seconds.
tp_im=zeros(max(obj),im_no);
for cel=1:max(obj)

 Ma=(Mask3==cel);
    mak1=gpuArray(Ma);
for ih=1:numbM

    if sum(sum(mak1(:,:,ih)))~=0 % imagesc(tp_im)  Elapsed time is 39.585907 seconds.
tp_im(cel,ih)=1;
    else
    end
   
end
cel
end
toc

%% split interrupted time series
tic
for cel=1:max(obj)
cel
tp_im2=diff(tp_im(cel,:));

tp1=find(tp_im2==1);
tp2=find(tp_im2==-1);

maxp=sum(sum(Mask3(:,:,numbM)==cel));


if length(tp1)==1 && length(tp2)==1  && maxp~=0% has one interruption

        for itx=tp1:numbM

        tp3=OAM_23121_tp3(Mask3(:,:,itx),cel,no_obj1,A);
       Mask3(:,:,itx)=tp3;

        end

    no_obj1=no_obj1+A;


elseif length(tp1)==1 && length(tp2)==1  && maxp==0% has one interruption

elseif length(tp1) == length(tp2)+1  && maxp~=0% 
       tp2(length(tp1))=numbM; 
     for itb=2:length(tp1) % starts at 2 because the first cell index remains unchanged 
        for itx=tp1(itb)+1:tp2(itb)
        
         tp3=OAM_23121_tp3(Mask3(:,:,itx),cel,no_obj1,A);
       Mask3(:,:,itx)=tp3;
        
        end
     no_obj1=no_obj1+A;
     end

elseif isempty(tp2) || isempty(tp1)  % its a normal cell, it's born and stays until the end

elseif length(tp1) == length(tp2) %!!! && sum(sum(masks{1,numbM}==cel))==0 % has one interruption !!! && sum(sum(masks{1,numbM}==cel))==0
     
if tp1(1)>tp2(1)
    tp2(length(tp2)+1)=numbM; 
     for itb=1:length(tp1) % starts at 2 because the first cell index remains unchanged 
        for itx=tp1(itb)+1:tp2(itb+1)
        
        tp3=OAM_23121_tp3(Mask3(:,:,itx),cel,no_obj1,A);
       Mask3(:,:,itx)=tp3;
        
        end
     no_obj1=no_obj1+A;
     end
elseif tp1(1)<tp2(1)

     for itb=2:length(tp1) % starts at 2 because the first cell index remains unchanged 
        for itx=tp1(itb)+1:tp2(itb)
        
      tp3=OAM_23121_tp3(Mask3(:,:,itx),cel,no_obj1,A);
       Mask3(:,:,itx)=tp3;
        
        end
     no_obj1=no_obj1+A;
     end

 elseif length(tp2)>1 % if it has multiple 

     for itb=2:length(tp1) % starts at 2 because the first cell index remains unchanged 
        for itx=tp1(itb)+1:tp2(itb)

        tp3=OAM_23121_tp3(Mask3(:,:,itx),cel,no_obj1,A);
       Mask3(:,:,itx)=tp3;

        end
     no_obj1=no_obj1+A;
     end

end

end

end
  toc

%% calculate size
no_obj2=max(unique(Mask3)); 
all_ob=zeros(no_obj2,im_no);
tic % 2.845733 seconds.
    for ccell=1:no_obj2 % ccell=1;
        ccell
 Maa=(Mask3==ccell);
    mak2=gpuArray(Maa);
 for io=1:im_no
pix=sum(sum(mak2(:,:,io)));
    all_ob(ccell,io)=pix;% figure;imagesc(all_ob)
 end
     end

  toc
%   figure(1); imagesc(all_ob);% title(num2str(cel))
    
    % cell existing birth and disappear times 
      cell_exists=zeros(2,size(all_ob,1));
    
    for itt2 = 1:size(all_ob,1)
    
        cell_exists(1:2,itt2)=[find(all_ob(itt2,:)~=0,1,'first'); find(all_ob(itt2,:)~=0,1,'last')];
   
    end

    no_obj=length(cell_exists);

    for its = 1:size(Mask3,3)
        Mask5{1,its} = Mask3(:,:,its);
    end
    clearvars Mask3

    Mask3 = Mask5;

save([sav_path pos '_ART_Track.mat'], 'all_ob','Mask3','no_obj','im_no','ccell2','cell_exists','shock_period','-v7.3');%,
