function diam_mat = diameter_of_matrix(diam_mat_model, N_bot_sat, N_samp, layer, conf_int, R_fiber, samp_pos_temp)
if diam_mat_model ~= 3
    diam_mat = {};%直径矩阵存储
else
    diam_mat = [];%直径矩阵存储
    diam_mat_sp = [];%临时、直观直径矩阵存储
end
if diam_mat_model == 1%模式一
    for i = 1:layer - 1
        N_bot_temp = N_bot_sat(i);
        diam_mat_temp = R_fiber + roundn(normrnd(0, conf_int/3, 1, N_bot_temp), -2);
        diam_mat_temp = repmat(diam_mat_temp, N_samp, 1);
        diam_mat{i, 1} = diam_mat_temp;
    end
else if diam_mat_model == 2%模式二
        for i = 1:layer - 1
            N_bot_temp = N_bot_sat(i);
            diam_mat_temp = R_fiber + roundn(normrnd(0, conf_int/3, N_samp, N_bot_temp), -2);
            diam_mat{i, 1} = diam_mat_temp;
        end
    else if diam_mat_model == 3%模式三，注：SP结尾代码特化
            for i = 1:length(samp_pos_temp)-1
                diam_mat_temp = R_fiber + roundn(normrnd(0, conf_int/3, 1, 1), -2);
                diam_mat_temp = repmat(diam_mat_temp, 1, N_samp);
                diam_mat_sp = [diam_mat_sp, diam_mat_temp];
            end
            
            for i = 1:length(samp_pos_temp)-1
                diam_mat_temp = (diam_mat_sp(N_samp*(i - 1)+1:N_samp*i))';
                diam_mat = [diam_mat, diam_mat_temp];
            end
        else if diam_mat_model == 4%模式四，注：SP结尾代码特化
                for i = 1:length(samp_pos_temp)-1
                    diam_mat_temp = R_fiber + roundn(normrnd(0, conf_int/3, 1, N_samp), -2);
                    diam_mat_sp = [diam_mat_sp, diam_mat_temp];
                end
                
                for i = 1:length(samp_pos_temp)-1
                    diam_mat_temp = (diam_mat_sp(N_samp*(i - 1)+1:N_samp*i))';
                    diam_mat = [diam_mat, diam_mat_temp];
                end
            end
        end
    end
end