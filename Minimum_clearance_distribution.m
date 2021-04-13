function ints_temp = Minimum_clearance_distribution(N_bot_temp, min_fiber, B_fiber, ints_mat_bot_rem, min_edge)
ints_temp = [];
ints_temp_all = 0;%已产生随即间隙的总长
ints_temp_t = 100;
nums = 1:N_bot_temp-1;%判断间隙中是否已产生随机数值
for i = 1:N_bot_temp-1
    while ints_temp_t > (min_fiber - B_fiber)||ints_temp_t < 0%产生的第一层的随机间隙矩阵的每个间隙距离大于某个定值
        ints_temp_t = roundn((ints_mat_bot_rem - ints_temp_all)*rand(1), -2);
    end
    kp = 0;
    while kp ~= 1
        pos = ceil((N_bot_temp - 1)*rand(1));
        pos_1 = (N_bot_temp - 1);
        if pos ~= 1||pos ~= pos_1
            kp = ismember(pos, nums);
            nums(pos) = NaN;
        end
    end
    ints_temp(pos, 1) = ints_temp_t + B_fiber;
    ints_temp_all = ints_temp_all + ints_temp_t;
    ints_temp_t = 100;
end
ints_all_bot_p = ints_mat_bot_rem - ints_temp_all;%判断上述代码随机后的剩余长度，直接进行均匀分配到每个间隙中
if ints_all_bot_p > 0
    ints_temp(1:end) = ints_temp(1:end) + ints_all_bot_p/(N_bot_temp - 1);
end
ints_temp = [min_edge, ints_temp', min_edge];