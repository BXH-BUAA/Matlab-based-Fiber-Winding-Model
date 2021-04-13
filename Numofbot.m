function N_bot = Numofbot(S_base,N_samp,R_fiber,B_fiber)
N_bot = 0;%底层匝数
fog_max = 0;%计算中光纤已占长度
while fog_max <= S_base
    N_bot = N_bot + 1;
    if N_bot>1
        fog_r_temp = R_fiber*ones(N_samp,N_bot);%底层光纤直径矩阵
        inter = B_fiber*ones(N_samp,N_bot-1);%底层光纤间隙矩阵
    else
        fog_r_temp = R_fiber*ones(N_samp,N_bot);
        inter = 0;
    end
    fog_sum = sum(fog_r_temp,2) + sum(inter,2);
    fog_max = max(fog_sum);
end
N_bot = N_bot-1;