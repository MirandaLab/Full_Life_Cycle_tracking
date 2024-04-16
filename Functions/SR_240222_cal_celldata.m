%% Function to calculate cell data - first occurance, last occurance, times the cell appears, number of times the cells disappears

function cell_data = cal_celldata(all_obj,ccel)
    cell_data = zeros(ccel,3);
    for iv=1:ccel
        cell_data(iv,1) = find(all_obj(iv,:) > 0,1, 'first'); %1st occurance
        cell_data(iv,2) = find(all_obj(iv,:) > 0,1, 'last'); %last occurance
    end
    for iv=1:ccel
        cell_data(iv,3) = cell_data(iv,2) - cell_data(iv,1) + 1; %times the cell appears !!!!!!! stress times
        aa1 = all_obj(iv,:);
        aa2 = aa1(cell_data(iv,1):cell_data(iv,2));
        aa3 = find(aa2 == 0);
        cell_data(iv,4) = numel(aa3); %number of times it disappears between 1st occurance and last occurance
        cell_data(iv,5) = (cell_data(iv,4)*100)/cell_data(iv,3); %percentage of times the cell disappers
    end
end