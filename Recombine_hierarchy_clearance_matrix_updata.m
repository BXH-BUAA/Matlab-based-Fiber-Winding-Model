function ints_mat_last = Recombine_hierarchy_clearance_matrix_updata(i, layer, ints_mat_last_temp,...
    T_base_left, T_base_left_ex, T_base_right, T_base_right_ex, N_bot_cal, diam_mat_last_temp)
if i ~= layer-1 || mod(i,2) == 0
    ints_mat_last = ints_mat_last_temp(:,(T_base_left + T_base_left_ex + 1 + 1):(T_base_left + T_base_left_ex + 1 + N_bot_cal));
    
    ints_mat_last_1 = sum(ints_mat_last_temp(:,1:T_base_left + T_base_left_ex + 1), 2) +...
        sum(diam_mat_last_temp(:,1:T_base_left + T_base_left_ex), 2);
    
    ints_mat_last_2 = sum(ints_mat_last_temp(:,((T_base_left + T_base_left_ex + 1) + N_bot_cal + 1):end), 2) +...
        sum(diam_mat_last_temp(:,((T_base_left + T_base_left_ex + N_bot_cal + 1) + 1):end), 2);
else if i == layer - 1&&mod(i,2) == 1
        ints_mat_last = ints_mat_last_temp(:,((end - (T_base_right + T_base_right_ex + 1)) - N_bot_cal + 1)...
            :(end - (T_base_right + T_base_right_ex + 1)));
        
        ints_mat_last_1 = sum(ints_mat_last_temp(:,1:((end - (T_base_right + T_base_right_ex + 1)) - N_bot_cal)), 2) +...
            sum(diam_mat_last_temp(:,1:((end - (T_base_right + T_base_right_ex) - N_bot_cal - 1))), 2);
        
        ints_mat_last_2 = sum(ints_mat_last_temp(:,((end - (T_base_right + T_base_right_ex + 1) + 1):end)), 2) +...
            sum(diam_mat_last_temp(:,(end-(T_base_right + T_base_right_ex) + 1):end), 2);
    end
end
ints_mat_last = [ints_mat_last_1, ints_mat_last, ints_mat_last_2];