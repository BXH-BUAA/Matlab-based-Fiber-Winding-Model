function fiber_core_bot = Bottom_core_matrix(diam_mat_bot, L_fiber_bot_max, ints_temp, N_bot_temp)
diam_mat_bot_max = diam_mat_bot(L_fiber_bot_max, :);
fiber_core_bot = [];%�洢�ײ���оλ�þ���
fiber_core_pos_last = ints_temp(1);
for i = 1:N_bot_temp
    fiber_core_pos = diam_mat_bot_max(i)/2 + fiber_core_pos_last;
    fiber_core_pos_last = fiber_core_pos + diam_mat_bot_max(i)/2 + ints_temp(i + 1);
    fiber_core_bot = [fiber_core_bot; fiber_core_pos];
end