function N_bot = Numofbot(S_base,N_samp,R_fiber,B_fiber)
N_bot = 0;%�ײ�����
fog_max = 0;%�����й�����ռ����
while fog_max <= S_base
    N_bot = N_bot + 1;
    if N_bot>1
        fog_r_temp = R_fiber*ones(N_samp,N_bot);%�ײ����ֱ������
        inter = B_fiber*ones(N_samp,N_bot-1);%�ײ���˼�϶����
    else
        fog_r_temp = R_fiber*ones(N_samp,N_bot);
        inter = 0;
    end
    fog_sum = sum(fog_r_temp,2) + sum(inter,2);
    fog_max = max(fog_sum);
end
N_bot = N_bot-1;