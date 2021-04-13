clc;
clear all;

%�ļ��´������Ƶ��ù��˻�Ϊ�ο����н�ģ
%ע�⣺��ģ�к��Ի��ѡ����ѵȶ�������ÿ�Ѽ�ΪԲ�����д����������ɹǼܸߴ����ʹ�

%�������йؼ�����������
%                   ������߲㼶��layer
%                   �㼶������N_bot_sat
%                   ��о����˻��Ǽܺ�����룺H_fiber_save
%                   �����ϲ���λ�ã�samp_pos
%                   ��϶����ints_mat
%                   ֱ������R_mat
%                   ��о����˻��Ǽܱ�����룺coc_mat
%                   �ײ���о����˻��Ǽ���ױ߾������fiber_core_bot

%����Ϊ��������
%ע�⣺�����г��ȡ�ֱ���Ȳ���Ӧ������um����;�ǶȵȲ���Ӧ�����ն�������

%���˻��Ǽܲ���
L_base = 5e5;%����
R_base = 1.2e5;%�ײ�ֱ��
A_base = 2;%��б�Ƕ�
S_base = L_base/cosd(A_base);%б�߳���
T_base_left = 1;%mÿ�����������ע�����ٷ�ʽΪ����1.5���ң�2.5
T_base_right = 2;

%���˲���
L_fiber = 25e9;%���˳���
R_fiber = 280;%g��������ֱ��
E_fiber = R_fiber*0.10;%���˼�ѹ�������ѹ�ٷֱ�
B_fiber = R_fiber*0.10;%���˻��ײ���С���
conf_int = R_fiber*10/1000;%����ֱ�����ƫ��
min_fiber = R_fiber*0.12;%����֮���϶��������
min_edge = 5;%�ײ���˾���˻��Ǽܵױ߾���

%��������
N_samp = 4;%ÿ��Բ���еĲ�������
D_samp = 360/N_samp;%�����Ƕȣ�����36������һȦ�а���10�����һ������

%�������λ�ü����һ�㼶��δ���ѣ����Ʋ���
samp_pos_control = 1;%���ڴ����м������λ�÷ǳ���ʱ������ͨ����������Ƿ�������λ�á�1���������λ�ã�0�����������λ��
N_bot_last_control = 0;%1���������һ�㼶��δ���ѣ���0�����������һ�㼶��δ���ѣ�
reconstructed_sg = 0;

%����������˻��ײ������������
N_bot_max_temp = Numofbot(S_base, N_samp, R_fiber, B_fiber);

%ע�⣺����N_bot_max_temp����������µ����ײ����������ڹ���ֱ��������������Ҫ����һ����Χ�ڵ��������
for N_bot = N_bot_max_temp-round(N_bot_max_temp/100):N_bot_max_temp-round(N_bot_max_temp/100)
    
    %ȷ���ײ������󣬺���ÿ�㼶����������֪�������԰��ոõײ������������м���
    
    %����������˻��Ĳ���λ��ͬʱ������˻���߲㼶
    %???�ص㻷��???
    %ע�⣺����������㼶���ڻᵼ����ʹ�õĹ��˳��ȳ���Ԥ��ֵ�����Լ�����Ĳ㼶��Ҫ��1��
    %      �����������һ�㼶��һ�����ѣ�������Ҫ�����������һ�㼶
    %      �����¼������λ�õĴ�����ʵ��ʹ���м����ʱ���ڲ���ʱ��Ӧ�ʵ��޸Ĳ�������������ʱ��
    %      ����ʱ��õ�ԭ��Ӧ����Ҫ���������ÿ����������ʼλ�õ�λ�ã��������������ࡢ���˳�ʱ���������󣬺�ʱ��
    S_base_rem = S_base - N_bot*R_fiber - (N_bot - 1)*B_fiber;%��ȥ����ֱ���Լ���׼��϶�⣬�ײ�ʣ�೤��
    ints_bot = B_fiber + S_base_rem/(N_bot - 1);%��������µĵײ��϶
    L_fiber_sat = 0;%ͳ����ʹ�ù��˳���
    L_fiber_sat_save = [];%�洢ÿ��ʹ�õĹ��˳���
    layer_temp = 1;%�㼶
    N_bot_sat = [];%�洢ÿ���㼶������
    H_fiber_save = {};%�洢���˽���Բ�ľ���Ǽܺ���߶�
    H_decr = sqrt(R_fiber^2 - ((R_fiber + ints_bot)/2)^2);%�����ڲ㼶֮��Ĺ��˺����Բ�ĸ߶Ȳ�
    while L_fiber_sat <= L_fiber
        H_fiber_layer = [];
        N_bot_idea = N_bot - (T_base_left + T_base_right + 1)*(layer_temp-1);%��ǰ�㼶�µ�����
        N_bot_sat = [N_bot_sat;N_bot_idea];
        if N_bot_idea <= 0
            break;
            disp('�����������?');
        end
        for i = 1:N_bot_idea
            H_fiber_layer_temp = R_base/2 - ((R_fiber/2)*sind(A_base) + (i - 1)*(R_fiber + ints_bot)*sind(A_base)) +...
                (R_fiber/2 + H_decr*(layer_temp - 1))*cosd(A_base) -...
                sind(A_base)*((layer_temp - 1)*(T_base_left*(R_fiber + ints_bot) + (R_fiber + ints_bot)/2));...
                %��ǰ�Ѿ���Ǽܺ���߶�
            H_fiber_layer = [H_fiber_layer, H_fiber_layer_temp];
        end
        H_fiber_save{layer_temp,1} = H_fiber_layer;
        
        %ͳ�Ƶ�ǰ�㼶���ù��˳���
        if mod(layer_temp,2)==1
            L_fiber_layer_sat = sum(2*pi*H_fiber_layer)+R_fiber*T_base_right;
        else
            L_fiber_layer_sat = sum(2*pi*H_fiber_layer)+R_fiber*T_base_left;
        end
        L_fiber_sat_save = [L_fiber_sat_save;L_fiber_layer_sat];
        L_fiber_sat = L_fiber_sat + L_fiber_layer_sat;
        
        layer_temp = layer_temp + 1;
        
    end
    N_bot_sat = N_bot_sat(1:end - 1);
    layer = layer_temp - 1;%���Ѳ㼶��Ϊlayer_temp - 2
    
    %���㵱ǰ�㼶�Ĺ��˾ݳ�ʼλ�õĲ�������
    samp_pos_temp = [0];%��������ͳ��
    samp_pos = [0];%������λ��ͳ��
    
    %����ǰlayer-1�����Ѳ㼶
    [samp_pos_temp, samp_pos] = Sampling_position_cal(samp_pos_control, layer, H_fiber_save,...
        samp_pos_temp, samp_pos, N_samp, reconstructed_sg);
    
    for mm=1:length(N_bot_sat)-1
        if mod(mm,2)==1
            samp_pos(sum(N_bot_sat(1:mm))*4+1:end)=samp_pos(sum(N_bot_sat(1:mm))*4:end)+R_fiber*T_base_right;
        else
            samp_pos(sum(N_bot_sat(1:mm))*4+1:end)=samp_pos(sum(N_bot_sat(1:mm))*4:end)+R_fiber*T_base_left;
        end
    end
    
    L_fiber_use = sum(L_fiber_sat_save(1:end-1));%���Ѳ㼶���ù��˳���
    
    if N_bot_last_control == 1
        %�������һ�㼶
        
        L_fiber_rme = L_fiber - L_fiber_use;%���һ�㼶���˳���
        %���һ�㼶��������
        H_fiber_layer_last = [];
        L_fiber_layer_last_sat = 0;
        i = 1;
        while L_fiber_layer_last_sat <= L_fiber_rme
            %�ж����һ�㼶����ż
            %ע�⣺�㼶����ż�����Ĺ��˴Ӹ����ͻ��Ǵӵ������ƣ�������Ҫ�ֱ��жϼ���
            if mod(layer,2) == 1
                %��ǰ�Ѿ���Ǽܺ���߶�
                H_fiber_layer_temp_last = R_base/2 - ((R_fiber/2)*sind(A_base) + (i - 1)*(R_fiber + ints_bot)*sind(A_base)) +...
                    (R_fiber/2 + H_decr*(layer - 1))*cosd(A_base) -...
                    sind(A_base)*((layer - 1)*(T_base_left*(R_fiber + ints_bot) + (R_fiber + ints_bot)/2));
            else
                %��ǰ�Ѿ���Ǽܺ���߶�
                H_fiber_layer_temp_last = R_base/2 - L_base*tand(A_base) +...
                    ((R_fiber/2)*sind(A_base) + (i - 1)*(R_fiber + ints_bot)*sind(A_base)) +...
                    (R_fiber/2 + H_decr*(layer - 1))*cosd(A_base) +...
                    sind(A_base)*((layer - 1)*(T_base_right*(R_fiber + ints_bot) + (R_fiber + ints_bot)/2));
            end
            H_fiber_layer_last = [H_fiber_layer_last, H_fiber_layer_temp_last];
            L_fiber_layer_last_sat = L_fiber_layer_last_sat + H_fiber_layer_temp_last*2*pi;
            i = i + 1;
        end
        N_bot_idea = i - 2;%���һ��δ���ѣ��ڸü����пɲ����Կ���
        N_bot_sat = [N_bot_sat;N_bot_idea];
        
        %������һ�㼶�Ĳ���λ��
        if samp_pos_control == 1
            if mod(layer,2) == 1
                S_fiber_layer_last = 2*pi*H_fiber_layer_last(1:end - 1)/N_samp;%ÿ�ѵĲ������
            else
                S_fiber_layer_last = 2*pi*fliplr(H_fiber_layer_last(1:end - 1))/N_samp;%ÿ�ѵĲ������
            end
            for m = 1:length(S_fiber_layer_last)
                for n = 1:N_samp
                    samp_pos = [samp_pos, samp_pos(end) + S_fiber_layer_last(m)];%��ò���λ��
                end
            end
        end
    end
    
    %�����ϴ��룬��֪����λ�ã���ɻ�ù���ֱ������
    %ע�⣺���ڸ�ֱ���������������������õ��ģ���ʵ������£�����λ���Ǵ���ƫ���
    %      ����ʵ���н�����ֱ����һ��һ�δ����ɼ���һ�����
    %      �ڸô����У�����û��ʵ�����ݣ���ͨ��������������ֱ������
    
    %����ֱ�����󣬷�Ϊ����ģʽ��һ��ÿ�㼶ÿ�ѹ���ֱ����ͬ��������ɣ�����ÿ���������ֱ����ȫ������ɣ�
    %                          ����ԭ��ͬһ����ȡ������һ�в�����λ�ã�ֱ���ԣ��ġ�ԭ��ͬ������
    diam_mat_model = 3;%ģʽ����
    diam_mat = diameter_of_matrix(diam_mat_model, N_bot_sat, N_samp, layer, conf_int, R_fiber, samp_pos_temp);%���ֱ������
    
    
    %���ֱ������֮����Ҫ��ͨ����ü�϶�������ӵײ���߲�������㣬�жϸ��㼶�Ĺ��˵�λ���Ƿ�����Ҫ��
    sat = 0;%�жϽ���Ƿ�����Ҫ��
    ints_mat = {};%��϶���󣬴洢ÿһ�㼶��϶��
    
    
    error_level = [];%�����б���㼶��¼
    ints_temp_save = [];
    while  sat == 0
        %������ʵ�ʲ��ƹ����У�ÿ��λ�õ�ϸ���ǲ���֪��Ҳ�޷���֪Ӧ�����������ĵײ��϶���ܻ���Ų����õĹ��˻�...
        %����Ŀǰ���ܹ�ͨ������ķ�ʽ������õײ�ļ�϶���󣬴Ӷ������ݸ߲�ǵļ�϶����
        %ע�⣺ʵ�����޷�ʹ�������ڹ��˵ļ�ౣ֤�ײ��϶�������Եײ�Ӧ���⿼�ǹ�����о��ļ��
        
        
        %ע�⣺���Ƶ�����ģʽ�£�����ÿ���㼶�Ĺ��˲���Ҫ����������˹Ǽܵĵױߣ����Կ��Կ����ڵײ��϶�����Լ��ײ���оλ�þ���...
        %      �����߸�����5um���ʵ����ȵ�����ֵ�������ļ��Գƹ����£����ڹ��˻��Ĳ��ֲ㼶����������˻��Ǽܵĵױߣ�����...
        %      �ײ��϶�����Լ��ײ���оλ�þ��������ֻ������Ϊ0
        
        
        %��õײ���оλ�þ����Լ��ײ��϶����
        N_bot_temp = N_bot_sat(1);
        
        %ע�⣺�ô�������������������Ļ�ȡ�ײ�ֱ������ʽ��ͬ������Ϊ�ô�����Ҫ�ص㿼�Ǳ��ֹ��˻���׶����̬��
        %      ���ܻ�������Ų㼶��������������һ������ǰ�����й涨�����ڲ㼶�仯T_base�������Ĺ涨��������Ҫͨ��������ʽ����õײ�ֱ������
        diam_mat_bot = diam_mat(:,1:N_bot_temp);
        
        %��ȥ���˻�ֱ���ͱ�Ҫ���˼���⣬ʣ�����С����
        ints_mat_bot_rem = S_base - max(sum(diam_mat_bot, 2)) - (N_bot_temp - 1)*B_fiber -10;%��10um�ǿ��ǵ�����
        [~,L_fiber_bot_max]= max(sum(diam_mat_bot, 2));
        
        %�������λ����ÿ����Сʣ�೤���µļ�϶����
        disp('�������λ����ÿ����Сʣ�೤���µļ�϶����');
        ints_temp = Minimum_clearance_distribution(N_bot_temp, min_fiber, B_fiber, ints_mat_bot_rem, min_edge);
        ints_temp_save = [ints_temp_save;ints_temp];
        
        %����ײ���оλ�þ���
        fiber_core_bot = Bottom_core_matrix(diam_mat_bot, L_fiber_bot_max, ints_temp, N_bot_temp);
        
        %����ײ��϶����
        ints_mat_bot = Bottom_clearance_matrix(N_samp, N_bot_temp, fiber_core_bot, diam_mat_bot, S_base);
        
        ints_mat{1,1} = ints_mat_bot;
        
        %��ʼ�Ʋ��ϲ㼶��϶����
        disp('��ʼ�Ʋ��ϲ㼶��϶����');
        
        delta_d_temp_layer_all=zeros(N_samp, N_bot);%�ײ���о��Ǽܱ�������ʼ����
        
        coc_mat = {};%��о��Ǽܱ������
        R_mat = {};%������Чֱ������洢
        N_bot_all = 0;
        reconstructed_signal = 0;
        for i = 1:layer-2
            ints_unroll_t = [];
            %ע�⣺�������һ�㼶������ż��δ֪��������Ҫ��������ۣ������һ�㼶Ϊ�棬����Ʒ����������ң��ɰ��յͲ㼶����
            N_bot_cal = N_bot_sat(i + 1);%��ǰ���㡢Ԥ�Ƶ���һ�㼶����
            N_bot_now = N_bot_sat(i);%��ǰ�㼶����
            N_bot_all = N_bot_all + N_bot_now;
            
            
            %�����������㼶֮����Լ���һ�������������Կ��Բ�����ǰ�����޵��ӵĹ���
            ints_mat_last_temp = ints_mat{i,1};
            diam_mat_last_temp = diam_mat(:,(N_bot_all - N_bot_now + 1):N_bot_all);
            
            if i == 1
                delta_d_temp_all=zeros(N_samp, N_bot_sat(i));
                R_mat{i, 1} = diam_mat_last_temp;
            end
            
            diam_mat_next = diam_mat(:,(N_bot_all + 1):N_bot_all + N_bot_cal);
            %���濪ʼ�жϵ�ǰ�������һ�㼶�ĵ�ǰ�����ĺ���ֱ��֮���Ƿ�����׶����״
            diam_mat_next_max = max(sum(diam_mat_next, 2) - N_bot_cal*R_fiber);
            diam_mat_next_max_t = diam_mat_next_max;
            k = 0;
            disp(diam_mat_next_max);
            while diam_mat_next_max_t > R_fiber/2
                k = k + 1;
                diam_mat_next_max_t = diam_mat_next_max - k*R_fiber;
            end
            T_base_left_ex = 0;
            T_base_right_ex = 0;
            if k ~= 0
                reconstructed_signal = 1;%��ʣ����˷������¹����źţ����¹�����������Լ�ֱ�����������
                if mod(i, 2) == 1
                    T_base_left_ex = k;
                    diam_mat_next = diam_mat_next(:,(k + 1):end);
                else
                    T_base_right_ex = k;
                    diam_mat_next = diam_mat_next(:,1:(end - k));
                end
                N_bot_sat((i + 1):end) = N_bot_sat((i + 1):end) - k;
                N_bot_cal = N_bot_cal - k;
            end
            R_mat{i + 1, 1} = diam_mat_next;
            
            %�ڴ˴���������������ۣ�
            %                     һ������߲㼶�������㼶
            %                     ������߲㼶Ϊ������
            %                     ������߲㼶Ϊż����
            %���У���һ��ڶ�������ڼ����϶����ʱ�����ڶ��Ǵ����Ƶ��ң����Կ���ͬ������
            if i ~= layer-1 || mod(i,2) == 0
                diam_mat_last = diam_mat_last_temp(:,(T_base_left + T_base_left_ex + 1):(T_base_left + T_base_left_ex + N_bot_cal + 1));
            else if i == layer - 1&&mod(i,2) == 1%�����������һ��δ���������㼶Ϊż���Ĺ��˲㼶
                    diam_mat_last = diam_mat_last_temp(:,((end - (T_base_right + T_base_right_ex) - N_bot_cal)...
                        :(end - (T_base_right + T_base_right_ex))));
                end
            end
            
            %����㼶��϶����
            ints_mat_last = Recombine_hierarchy_clearance_matrix_updata(i, layer,...
                ints_mat_last_temp, T_base_left, T_base_left_ex, T_base_right, T_base_right_ex, N_bot_cal, diam_mat_last_temp);
            
            %����ÿ�㼶��϶����ÿ�зֱ���д���
            for m = 1:N_samp
                if i ~= layer-1 || mod(i,2) == 0
                    delta_d_temp = (delta_d_temp_all(m,(T_base_left + T_base_left_ex + 1)...
                        :(T_base_left + T_base_left_ex + N_bot_cal + 1)))';
                else if i == layer - 1&&mod(i,2) == 1
                        delta_d_temp = (delta_d_temp_all(m,((end - T_base_right + T_base_right_ex - N_bot_cal)...
                            :(end - T_base_right + T_base_right_ex))))';
                    end
                end
                ints_temp = ints_mat_last(m,:);%��ȡ�Ͳ㼶m�м�϶����
                data_unroll_temp_1 = diam_mat_last(m,:);%��ȡ�Ͳ㼶m��ֱ������
                data_unroll_temp_2 = diam_mat_next(m,:);%��ȡ�߲㼶m��ֱ������
                ints_unroll_temp = [];
                ints_unroll_temp_2 = 0;
                
                %ע�⣺ÿ��������һ�μ�϶�����Ͳ㼶���ˣ�A\B����һ�߲㼶���ˣ�C����ɣ���о���Ӻ�Ϊ�����Σ����������ι�ʽ���м���
                
                for j = 1:N_bot_cal%ÿ��ÿ�ѵ��Ʋ���һ���ĵ�����϶���
                    %��ȡ��Ŀ���������˾���˻��Ǽܵײ�����
                    if j ~= 1
                        ints_1 = sum(ints_temp(1:j)) + sum(data_unroll_temp_1(1:j - 1));
                    else
                        ints_1 = sum(ints_temp(1:j));
                    end
                    
                    %��ȡ��Ŀ��������������˼��
                    ints_2 = ints_temp(j+1);
                    
                    %��ȡ��Ŀ�����
                    fog_r_unroll_1 = data_unroll_temp_1(j);
                    fog_r_unroll_2 = data_unroll_temp_1(j+1);
                    fog_r_unroll_3 = data_unroll_temp_2(j);
                    
                    a = (fog_r_unroll_3 + fog_r_unroll_1)/2;%A\C������оֱ�߾���
                    %dΪA\B�����˵���о�߶Ȳ�
                    if i == 1
                        d = (delta_d_temp(j,1) + fog_r_unroll_1 - delta_d_temp(j + 1,1) - fog_r_unroll_2)/2;
                    else
                        d = delta_d_temp(j,1) - delta_d_temp(j + 1,1);
                    end
                    b = sqrt(((fog_r_unroll_2 + fog_r_unroll_1)/2 + ints_2)^2 + d^2);%A\B��������о����������
                    c = (fog_r_unroll_3 + fog_r_unroll_2)/2;%B\C������оֱ�߾���
                    theta1 = rad2deg(acos((a^2 + b^2 - c^2)/(2*a*b)));%A\C��A\B�������߹��ɵĽǶ�
                    theta2 = rad2deg(asin(d/b));%A\B���߹��ɵ������ĽǶ�
                    theta = theta1-theta2;%A\C���߹��ɵ������ĽǶ�
                    delta_d = a*cosd(theta);%A\C��������о���߾����ں����ϵķ���
                    if i == 1
                        delta_d_temp(j,1) = delta_d_temp(j,1) + a*sind(theta) + fog_r_unroll_1/2;
                    else
                        delta_d_temp(j,1) = delta_d_temp(j,1) + a*sind(theta);%�߲㼶��о��Ǽܱ������
                    end
                    %�����϶
                    ints_unroll_temp_1 = ints_1+fog_r_unroll_1/2 + delta_d - fog_r_unroll_3/2 - ints_unroll_temp_2;
                    if ints_unroll_temp_1 < -1*E_fiber
                        break%�жϸü�϶����Ƿ������϶�жϻ�׼��������ֱ�ӽ��������¼���
                        disp('��϶δ��������СҪ��');
                    end
                    ints_unroll_temp = [ints_unroll_temp;ints_unroll_temp_1];%�洢�߲㼶��϶����
                    %��ȡ�߲㼶�����ұ߾�Ǽܵױ߾���
                    ints_unroll_temp_2 = sum(ints_unroll_temp) + sum(data_unroll_temp_2(1:j));
                end
                if ints_unroll_temp_1 < -1*E_fiber
                    break%����
                end
                delta_d_temp_all(m,1:N_bot_cal) = (delta_d_temp(1:N_bot_cal))';
                ints_unroll_temp_1_last = S_base - sum(ints_unroll_temp) - sum(data_unroll_temp_2);%���������˾�Ǽܵײ�����
                ints_unroll_temp = [ints_unroll_temp;ints_unroll_temp_1_last];
                ints_unroll_t = [ints_unroll_t;ints_unroll_temp'];%�����׼�󣬽����м�϶���ݴ洢
            end
            if ints_unroll_temp_1 < -1*E_fiber
                disp(['����㼶��', num2str(i)]);
                disp(['�����϶��', num2str(ints_unroll_temp_1)]);
                disp(['����������', num2str(j)]);
                error_level=[error_level,i];
                %tabulate(error_level);
                break%����
            end
            ints_mat{i + 1,1} = ints_unroll_t;%�洢��϶����
            if i == 1
                coc_mat{i,1} = (R_mat{i,1})/2;
            end
            delta_d_temp_all = delta_d_temp_all(:,1:N_bot_cal);
            coc_mat{i + 1,1}=delta_d_temp_all;%�洢��о��Ǽܱ������
            
            if reconstructed_signal == 1
                disp('Data Reconstructed ?');
                turns = sum(N_bot_sat(1:(i + 1)));%�����ع��źź����ù��˳���
                samp_pos = samp_pos(1:turns*N_samp+1);
                samp_pos_temp = samp_pos_temp(1:turns+1);
                reconstructed_sg = i + 1;
                [diam_mat, samp_pos, samp_pos_temp] = reconstructed_fiber_test(i, samp_pos_control, layer, H_fiber_save, R_fiber, A_base, ints_bot, R_base, H_decr,...
                    samp_pos_temp, samp_pos, N_samp, reconstructed_sg,...
                    N_bot_sat, T_base_left, T_base_left_ex, T_base_right, conf_int, L_fiber_sat_save, L_fiber);
            end
            
        end
        if i == layer - 2
            sat = 1;%���Ʋ������в㼶�󣬸���sat���в�����ֹ�Ʋ�
        end
    end
end
