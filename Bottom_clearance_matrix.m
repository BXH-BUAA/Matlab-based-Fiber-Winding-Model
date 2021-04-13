function ints_mat_bot = Bottom_clearance_matrix(N_samp, N_bot_temp, fiber_core_bot, diam_mat_bot, S_base)
ints_mat_bot = [];
fiber_core_bot_temp = 0;
diam_mat_bot_temp = zeros(N_samp, 1);
for j = 1:N_bot_temp
    %��ȡ�ײ��϶��ÿ�����ݵļ�϶����
    ints_mat_bot(:, j) = fiber_core_bot(j) - diam_mat_bot(:, j)/2 - fiber_core_bot_temp - diam_mat_bot_temp;
    fiber_core_bot_temp = fiber_core_bot(j);
    diam_mat_bot_temp = diam_mat_bot(:, j)/2;
end
%��ȡ�ײ��϶�������һ������˻��Ǽ�֮��ļ�϶
ints_mat_bot(:, j+1) = S_base - fiber_core_bot_temp - diam_mat_bot_temp;