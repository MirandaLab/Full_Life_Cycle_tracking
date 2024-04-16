
function [IS1B]=OAM_231216_bina(IS1)
IS1B = IS1;% figure;imagesc(IS1); title('IS1');
IS1B(IS1~=0)=1;% figure;imagesc(IS1B) - mask of the tracked past in 1s and 0s