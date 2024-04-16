%% Function to calculate sizes of all the cells tracked in the time series

function all_obj = cal_allob(ccel,TETC,rang)

    all_obj = zeros(ccel,size(TETC,2));
    for iv=1:ccel
        for its=rang
            if ~isempty(TETC{1,its})
                all_obj(iv,its) = sum(sum(logical(TETC{1,its} == iv)));
            else
                all_obj(iv,its) = -1;
            end
        end
    end


end