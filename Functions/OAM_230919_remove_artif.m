function [I] = OAM_230919_remove_artif(I2A)


    I2B = imopen(logical(I2A), strel('disk', 6)); % 3 for 752x1024 and 6 for 1024x2048
    I2C = imdilate(I2B, strel('disk', 6)); % 3 for 752x1024 and 6 for 1024x2048

    I3 = uint16(I2A).*uint16(I2C);
    objs = unique(I3);

    I4 = zeros(size(I2A));

    for i = 1:size(objs,1)
        I4(I2A == objs(i)) = objs(i);
    end
    % I4 is the new artifacts free image
    
    I=uint16(I4);


% 