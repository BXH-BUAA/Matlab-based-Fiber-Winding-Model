function ints_temp = Minimum_clearance_distribution(N_bot_temp, min_fiber, B_fiber, ints_mat_bot_rem, min_edge)
ints_temp = [];
ints_temp_all = 0;%�Ѳ����漴��϶���ܳ�
ints_temp_t = 100;
nums = 1:N_bot_temp-1;%�жϼ�϶���Ƿ��Ѳ��������ֵ
for i = 1:N_bot_temp-1
    while ints_temp_t > (min_fiber - B_fiber)||ints_temp_t < 0%�����ĵ�һ��������϶�����ÿ����϶�������ĳ����ֵ
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
ints_all_bot_p = ints_mat_bot_rem - ints_temp_all;%�ж���������������ʣ�೤�ȣ�ֱ�ӽ��о��ȷ��䵽ÿ����϶��
if ints_all_bot_p > 0
    ints_temp(1:end) = ints_temp(1:end) + ints_all_bot_p/(N_bot_temp - 1);
end
ints_temp = [min_edge, ints_temp', min_edge];