function [samp_pos_tmep, samp_pos] = Sampling_position_cal(samp_pos_control, layer, H_fiber_save, samp_pos_tmep, samp_pos, N_samp, reconstructed_sg)
if samp_pos_control == 1
    for layer_temp = 1 + reconstructed_sg:layer-1
        H_fiber_layer_pos = H_fiber_save{layer_temp,1};
        if mod(layer_temp,2) == 0
            S_fiber_layer = fliplr(H_fiber_layer_pos);
        else
            S_fiber_layer = H_fiber_layer_pos;
        end
        S_fiber_layer = 2*pi*S_fiber_layer;%每匝的长度
        for m = 1:length(S_fiber_layer)
            samp_pos_tmep = [samp_pos_tmep, samp_pos_tmep(end) + S_fiber_layer(m)];%获得每匝采样长度
        end
    end
end
%下面将每匝采样进行具体细分处理
if N_samp == 4
    for m = 1:length(samp_pos_tmep)-1
        K4 = samp_pos_tmep(m + 1);
        K1 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*1/N_samp;
        K2 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*2/N_samp;
        K3 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*3/N_samp;
        samp_pos(1 + (m*N_samp - 3)) = K1;
        samp_pos(1 + (m*N_samp - 2)) = K2;
        samp_pos(1 + (m*N_samp - 1)) = K3;
        samp_pos(1 + (m*N_samp - 0)) = K4;
    end
else if N_samp == 5
        for m = 1:length(samp_pos_tmep)-1
            K5 = samp_pos_tmep(m + 1);
            K1 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*1/N_samp;
            K2 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*2/N_samp;
            K3 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*3/N_samp;
            K4 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*4/N_samp;
                samp_pos(1 + (m*N_samp - 4)) = K1;
                samp_pos(1 + (m*N_samp - 3)) = K2;
                samp_pos(1 + (m*N_samp - 2)) = K3;
                samp_pos(1 + (m*N_samp - 1)) = K4;
                samp_pos(1 + (m*N_samp - 0)) = K5;
        end
    else if N_samp == 6
            for m = 1:length(samp_pos_tmep)-1
                K6 = samp_pos_tmep(m + 1);
                K1 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*1/N_samp;
                K2 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*2/N_samp;
                K3 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*3/N_samp;
                K4 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*4/N_samp;
                K5 = samp_pos_tmep(m) + (samp_pos_tmep(m + 1) - samp_pos_tmep(m))*5/N_samp;
                samp_pos(1 + (m*N_samp - 5)) = K1;
                samp_pos(1 + (m*N_samp - 4)) = K2;
                samp_pos(1 + (m*N_samp - 3)) = K3;
                samp_pos(1 + (m*N_samp - 2)) = K4;
                samp_pos(1 + (m*N_samp - 1)) = K5;
                samp_pos(1 + (m*N_samp - 0)) = K6;
            end
        end
    end
end