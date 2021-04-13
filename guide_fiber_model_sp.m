clc;
clear all;

%文件下代码以制导用光纤环为参考进行建模
%注意：建模中忽略换匝、跨匝等动作，将每匝简化为圆环进行处理，处理方向由骨架高处往低处

%本代码中关键变量包括：
%                   理论最高层级：layer
%                   层级匝数：N_bot_sat
%                   纤芯距光纤环骨架横轴距离：H_fiber_save
%                   光纤上采样位置：samp_pos
%                   间隙矩阵：ints_mat
%                   直径矩阵：R_mat
%                   纤芯距光纤环骨架表面距离：coc_mat
%                   底层纤芯距光纤环骨架左底边距离矩阵：fiber_core_bot

%以下为基本参数
%注意：代码中长度、直径等参数应均按照um处理;角度等参数应均按照度数处理

%光纤环骨架参数
L_base = 5e5;%长度
R_base = 1.2e5;%底部直径
A_base = 2;%倾斜角度
S_base = L_base/cosd(A_base);%斜边长度
T_base_left = 1;%m每层减少匝数，注：减少方式为：左：1.5，右：2.5
T_base_right = 2;

%光纤参数
L_fiber = 25e9;%光纤长度
R_fiber = 280;%g光纤理想直径
E_fiber = R_fiber*0.10;%光纤挤压最大允许挤压百分比
B_fiber = R_fiber*0.10;%光纤环底层最小间距
conf_int = R_fiber*10/1000;%光纤直径最大偏差
min_fiber = R_fiber*0.12;%光纤之间间隙初步控制
min_edge = 5;%底层光纤距光纤环骨架底边距离

%采样参数
N_samp = 4;%每个圆环中的采样个数
D_samp = 360/N_samp;%采样角度，例：36代表在一圈中按照10°进行一个采样

%计算采样位置及最后一层级（未满匝）控制参数
samp_pos_control = 1;%由于代码中计算采样位置非常耗时，于是通过代码控制是否计算采样位置。1即计算采样位置，0即不计算采样位置
N_bot_last_control = 0;%1即计算最后一层级（未满匝），0即不计算最后一层级（未满匝）
reconstructed_sg = 0;

%计算理想光纤环底层最大允许匝数
N_bot_max_temp = Numofbot(S_base, N_samp, R_fiber, B_fiber);

%注意：由于N_bot_max_temp是理想情况下的最大底层匝数，由于光纤直径存在误差，所以需要考虑一定范围内的匝数误差
for N_bot = N_bot_max_temp-round(N_bot_max_temp/100):N_bot_max_temp-round(N_bot_max_temp/100)
    
    %确定底层匝数后，后续每层级匝数都是已知数，所以按照该底层匝数继续进行计算
    
    %计算理想光纤环的采样位置同时计算光纤环最高层级
    %???重点环节???
    %注意：计算出的最后层级由于会导致所使用的光纤长度超过预定值，所以计算出的层级需要减1，
    %      并且由于最后一层级不一定满匝，所以需要单独考虑最后一层级
    %      另：以下计算采样位置的代码在实际使用中极其耗时，在测试时，应适当修改参数，减少运行时间
    %      运行时间久的原因应是需要计算光纤上每个采样点距初始位置的位置，当采样点数量多、光纤长时，运算量大，耗时久
    S_base_rem = S_base - N_bot*R_fiber - (N_bot - 1)*B_fiber;%除去光纤直径以及标准间隙外，底层剩余长度
    ints_bot = B_fiber + S_base_rem/(N_bot - 1);%理想情况下的底层间隙
    L_fiber_sat = 0;%统计已使用光纤长度
    L_fiber_sat_save = [];%存储每层使用的光纤长度
    layer_temp = 1;%层级
    N_bot_sat = [];%存储每个层级的匝数
    H_fiber_save = {};%存储光纤截面圆心距离骨架横轴高度
    H_decr = sqrt(R_fiber^2 - ((R_fiber + ints_bot)/2)^2);%两相邻层级之间的光纤横截面圆心高度差
    while L_fiber_sat <= L_fiber
        H_fiber_layer = [];
        N_bot_idea = N_bot - (T_base_left + T_base_right + 1)*(layer_temp-1);%当前层级下的匝数
        N_bot_sat = [N_bot_sat;N_bot_idea];
        if N_bot_idea <= 0
            break;
            disp('参数输入错误！?');
        end
        for i = 1:N_bot_idea
            H_fiber_layer_temp = R_base/2 - ((R_fiber/2)*sind(A_base) + (i - 1)*(R_fiber + ints_bot)*sind(A_base)) +...
                (R_fiber/2 + H_decr*(layer_temp - 1))*cosd(A_base) -...
                sind(A_base)*((layer_temp - 1)*(T_base_left*(R_fiber + ints_bot) + (R_fiber + ints_bot)/2));...
                %当前匝距离骨架横轴高度
            H_fiber_layer = [H_fiber_layer, H_fiber_layer_temp];
        end
        H_fiber_save{layer_temp,1} = H_fiber_layer;
        
        %统计当前层级所用光纤长度
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
    layer = layer_temp - 1;%满匝层级数为layer_temp - 2
    
    %计算当前层级的光纤据初始位置的采样距离
    samp_pos_temp = [0];%采样点间距统计
    samp_pos = [0];%采样点位置统计
    
    %考虑前layer-1个满匝层级
    [samp_pos_temp, samp_pos] = Sampling_position_cal(samp_pos_control, layer, H_fiber_save,...
        samp_pos_temp, samp_pos, N_samp, reconstructed_sg);
    
    for mm=1:length(N_bot_sat)-1
        if mod(mm,2)==1
            samp_pos(sum(N_bot_sat(1:mm))*4+1:end)=samp_pos(sum(N_bot_sat(1:mm))*4:end)+R_fiber*T_base_right;
        else
            samp_pos(sum(N_bot_sat(1:mm))*4+1:end)=samp_pos(sum(N_bot_sat(1:mm))*4:end)+R_fiber*T_base_left;
        end
    end
    
    L_fiber_use = sum(L_fiber_sat_save(1:end-1));%满匝层级所用光纤长度
    
    if N_bot_last_control == 1
        %考虑最后一层级
        
        L_fiber_rme = L_fiber - L_fiber_use;%最后一层级光纤长度
        %最后一层级可绕匝数
        H_fiber_layer_last = [];
        L_fiber_layer_last_sat = 0;
        i = 1;
        while L_fiber_layer_last_sat <= L_fiber_rme
            %判断最后一层级数奇偶
            %注意：层级的奇偶决定的光纤从高往低还是从低往高绕，所以需要分别判断计算
            if mod(layer,2) == 1
                %当前匝距离骨架横轴高度
                H_fiber_layer_temp_last = R_base/2 - ((R_fiber/2)*sind(A_base) + (i - 1)*(R_fiber + ints_bot)*sind(A_base)) +...
                    (R_fiber/2 + H_decr*(layer - 1))*cosd(A_base) -...
                    sind(A_base)*((layer - 1)*(T_base_left*(R_fiber + ints_bot) + (R_fiber + ints_bot)/2));
            else
                %当前匝距离骨架横轴高度
                H_fiber_layer_temp_last = R_base/2 - L_base*tand(A_base) +...
                    ((R_fiber/2)*sind(A_base) + (i - 1)*(R_fiber + ints_bot)*sind(A_base)) +...
                    (R_fiber/2 + H_decr*(layer - 1))*cosd(A_base) +...
                    sind(A_base)*((layer - 1)*(T_base_right*(R_fiber + ints_bot) + (R_fiber + ints_bot)/2));
            end
            H_fiber_layer_last = [H_fiber_layer_last, H_fiber_layer_temp_last];
            L_fiber_layer_last_sat = L_fiber_layer_last_sat + H_fiber_layer_temp_last*2*pi;
            i = i + 1;
        end
        N_bot_idea = i - 2;%最后一匝未满匝，在该计算中可不予以考虑
        N_bot_sat = [N_bot_sat;N_bot_idea];
        
        %获得最后一层级的采样位置
        if samp_pos_control == 1
            if mod(layer,2) == 1
                S_fiber_layer_last = 2*pi*H_fiber_layer_last(1:end - 1)/N_samp;%每匝的采样间距
            else
                S_fiber_layer_last = 2*pi*fliplr(H_fiber_layer_last(1:end - 1))/N_samp;%每匝的采样间距
            end
            for m = 1:length(S_fiber_layer_last)
                for n = 1:N_samp
                    samp_pos = [samp_pos, samp_pos(end) + S_fiber_layer_last(m)];%获得采样位置
                end
            end
        end
    end
    
    %由以上代码，已知采样位置，则可获得光纤直径矩阵
    %注意：由于该直径矩阵是由理想情况计算得到的，而实际情况下，采样位置是存在偏差的
    %      所以实际中将光纤直径按一段一段处理，可减少一定误差
    %      在该代码中，由于没有实际数据，将通过随机生成来获得直径矩阵
    
    %生成直径矩阵，分为两种模式，一、每层级每匝光纤直径相同，随机生成；二、每个采样点的直径完全随机生成；
    %                          三、原理同一，获取光纤上一行采样点位置，直观性；四、原理同二及三
    diam_mat_model = 3;%模式控制
    diam_mat = diameter_of_matrix(diam_mat_model, N_bot_sat, N_samp, layer, conf_int, R_fiber, samp_pos_temp);%获得直径矩阵
    
    
    %获得直径矩阵之后，需要再通过获得间隙矩阵来从底层向高层进行推算，判断各层级的光纤的位置是否满足要求
    sat = 0;%判断结果是否满足要求
    ints_mat = {};%间隙矩阵，存储每一层级间隙用
    
    
    error_level = [];%推算中报错层级记录
    ints_temp_save = [];
    while  sat == 0
        %由于在实际缠绕过程中，每个位置的细节是不可知，也无法获知应当适用怎样的底层间隙才能获得排布良好的光纤环...
        %所以目前仅能够通过随机的方式，来获得底层的间隙矩阵，从而来推演高层记的间隙矩阵
        %注意：实际中无法使两根相邻光纤的间距保证底层间隙矩阵，所以底层应额外考虑光纤纤芯间的间距
        
        
        %注意：在制导光纤模式下，由于每个层级的光纤不需要必须紧贴光纤骨架的底边，所以可以考虑在底层间隙矩阵以及底层纤芯位置矩阵...
        %      的两边各增加5um或适当长度的冗余值；而在四级对称光纤下，由于光纤环的部分层级必须紧贴光纤环骨架的底边，所以...
        %      底层间隙矩阵以及底层纤芯位置矩阵的两边只能设置为0
        
        
        %获得底层纤芯位置矩阵以及底层间隙矩阵
        N_bot_temp = N_bot_sat(1);
        
        %注意：该代码中与其余两个代码的获取底层直径矩阵方式不同，是因为该代码需要重点考虑保持光纤环的锥形形态，
        %      可能会存在随着层级的上升，匝数不一定按照前代码中规定的相邻层级变化T_base的匝数的规定，所以需要通过下述方式来获得底层直径矩阵
        diam_mat_bot = diam_mat(:,1:N_bot_temp);
        
        %除去光纤环直径和必要光纤间距外，剩余的最小长度
        ints_mat_bot_rem = S_base - max(sum(diam_mat_bot, 2)) - (N_bot_temp - 1)*B_fiber -10;%减10um是考虑到冗余
        [~,L_fiber_bot_max]= max(sum(diam_mat_bot, 2));
        
        %计算采样位置下每行最小剩余长度下的间隙分配
        disp('计算采样位置下每行最小剩余长度下的间隙分配');
        ints_temp = Minimum_clearance_distribution(N_bot_temp, min_fiber, B_fiber, ints_mat_bot_rem, min_edge);
        ints_temp_save = [ints_temp_save;ints_temp];
        
        %计算底层纤芯位置矩阵
        fiber_core_bot = Bottom_core_matrix(diam_mat_bot, L_fiber_bot_max, ints_temp, N_bot_temp);
        
        %计算底层间隙矩阵
        ints_mat_bot = Bottom_clearance_matrix(N_samp, N_bot_temp, fiber_core_bot, diam_mat_bot, S_base);
        
        ints_mat{1,1} = ints_mat_bot;
        
        %开始推测上层级间隙矩阵
        disp('开始推测上层级间隙矩阵');
        
        delta_d_temp_layer_all=zeros(N_samp, N_bot);%底层纤芯距骨架表面距离初始矩阵
        
        coc_mat = {};%纤芯距骨架表面距离
        R_mat = {};%光纤有效直径矩阵存储
        N_bot_all = 0;
        reconstructed_signal = 0;
        for i = 1:layer-2
            ints_unroll_t = [];
            %注意：由于最后一层级级数奇偶性未知，所以需要分情况讨论，若最后一层级为奇，其缠绕方向由左至右，可按照低层级处理
            N_bot_cal = N_bot_sat(i + 1);%当前计算、预计的下一层级匝数
            N_bot_now = N_bot_sat(i);%当前层级匝数
            N_bot_all = N_bot_all + N_bot_now;
            
            
            %由于相邻两层级之间相对减少一定量匝数，所以可以不考虑前后几匝无叠加的光纤
            ints_mat_last_temp = ints_mat{i,1};
            diam_mat_last_temp = diam_mat(:,(N_bot_all - N_bot_now + 1):N_bot_all);
            
            if i == 1
                delta_d_temp_all=zeros(N_samp, N_bot_sat(i));
                R_mat{i, 1} = diam_mat_last_temp;
            end
            
            diam_mat_next = diam_mat(:,(N_bot_all + 1):N_bot_all + N_bot_cal);
            %下面开始判断当前计算的下一层级的当前匝数的横向直径之和是否满足锥形形状
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
                reconstructed_signal = 1;%对剩余光纤发出重新构建信号，重新构建采样间距以及直径矩阵等数据
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
            
            %在此处共分三种情况讨论：
            %                     一、除最高层级外的其余层级
            %                     二、最高层级为奇数层
            %                     三、最高层级为偶数层
            %其中，第一与第二中情况在计算间隙矩阵时，由于都是从左绕到右，所以可以同类讨论
            if i ~= layer-1 || mod(i,2) == 0
                diam_mat_last = diam_mat_last_temp(:,(T_base_left + T_base_left_ex + 1):(T_base_left + T_base_left_ex + N_bot_cal + 1));
            else if i == layer - 1&&mod(i,2) == 1%单独考虑最后一层未满的且最后层级为偶数的光纤层级
                    diam_mat_last = diam_mat_last_temp(:,((end - (T_base_right + T_base_right_ex) - N_bot_cal)...
                        :(end - (T_base_right + T_base_right_ex))));
                end
            end
            
            %重组层级间隙矩阵
            ints_mat_last = Recombine_hierarchy_clearance_matrix_updata(i, layer,...
                ints_mat_last_temp, T_base_left, T_base_left_ex, T_base_right, T_base_right_ex, N_bot_cal, diam_mat_last_temp);
            
            %按照每层级间隙矩阵每行分别进行处理
            for m = 1:N_samp
                if i ~= layer-1 || mod(i,2) == 0
                    delta_d_temp = (delta_d_temp_all(m,(T_base_left + T_base_left_ex + 1)...
                        :(T_base_left + T_base_left_ex + N_bot_cal + 1)))';
                else if i == layer - 1&&mod(i,2) == 1
                        delta_d_temp = (delta_d_temp_all(m,((end - T_base_right + T_base_right_ex - N_bot_cal)...
                            :(end - T_base_right + T_base_right_ex))))';
                    end
                end
                ints_temp = ints_mat_last(m,:);%获取低层级m行间隙矩阵
                data_unroll_temp_1 = diam_mat_last(m,:);%获取低层级m行直径矩阵
                data_unroll_temp_2 = diam_mat_next(m,:);%获取高层级m行直径矩阵
                ints_unroll_temp = [];
                ints_unroll_temp_2 = 0;
                
                %注意：每次推算下一次间隙由两低层级光纤（A\B）和一高层级光纤（C）组成，纤芯连接后为三角形，则由三角形公式进行计算
                
                for j = 1:N_bot_cal%每匝每匝的推测下一级的单个间隙宽度
                    %获取三目标光纤左光纤距光纤环骨架底部距离
                    if j ~= 1
                        ints_1 = sum(ints_temp(1:j)) + sum(data_unroll_temp_1(1:j - 1));
                    else
                        ints_1 = sum(ints_temp(1:j));
                    end
                    
                    %获取三目标光纤左右两光纤间距
                    ints_2 = ints_temp(j+1);
                    
                    %获取三目标光纤
                    fog_r_unroll_1 = data_unroll_temp_1(j);
                    fog_r_unroll_2 = data_unroll_temp_1(j+1);
                    fog_r_unroll_3 = data_unroll_temp_2(j);
                    
                    a = (fog_r_unroll_3 + fog_r_unroll_1)/2;%A\C光纤纤芯直线距离
                    %d为A\B两光纤的纤芯高度差
                    if i == 1
                        d = (delta_d_temp(j,1) + fog_r_unroll_1 - delta_d_temp(j + 1,1) - fog_r_unroll_2)/2;
                    else
                        d = delta_d_temp(j,1) - delta_d_temp(j + 1,1);
                    end
                    b = sqrt(((fog_r_unroll_2 + fog_r_unroll_1)/2 + ints_2)^2 + d^2);%A\B两光纤纤芯距离横向分量
                    c = (fog_r_unroll_3 + fog_r_unroll_2)/2;%B\C光纤纤芯直线距离
                    theta1 = rad2deg(acos((a^2 + b^2 - c^2)/(2*a*b)));%A\C、A\B两条连线构成的角度
                    theta2 = rad2deg(asin(d/b));%A\B连线构成的与横轴的角度
                    theta = theta1-theta2;%A\C连线构成的与横轴的角度
                    delta_d = a*cosd(theta);%A\C两光纤纤芯连线距离在横向上的分量
                    if i == 1
                        delta_d_temp(j,1) = delta_d_temp(j,1) + a*sind(theta) + fog_r_unroll_1/2;
                    else
                        delta_d_temp(j,1) = delta_d_temp(j,1) + a*sind(theta);%高层级纤芯距骨架表面距离
                    end
                    %计算间隙
                    ints_unroll_temp_1 = ints_1+fog_r_unroll_1/2 + delta_d - fog_r_unroll_3/2 - ints_unroll_temp_2;
                    if ints_unroll_temp_1 < -1*E_fiber
                        break%判断该间隙宽度是否满足间隙判断基准，不满足直接结束，重新计算
                        disp('间隙未能满足最小要求！');
                    end
                    ints_unroll_temp = [ints_unroll_temp;ints_unroll_temp_1];%存储高层级间隙矩阵
                    %获取高层级光纤右边距骨架底边距离
                    ints_unroll_temp_2 = sum(ints_unroll_temp) + sum(data_unroll_temp_2(1:j));
                end
                if ints_unroll_temp_1 < -1*E_fiber
                    break%跳出
                end
                delta_d_temp_all(m,1:N_bot_cal) = (delta_d_temp(1:N_bot_cal))';
                ints_unroll_temp_1_last = S_base - sum(ints_unroll_temp) - sum(data_unroll_temp_2);%计算最后光纤距骨架底部距离
                ints_unroll_temp = [ints_unroll_temp;ints_unroll_temp_1_last];
                ints_unroll_t = [ints_unroll_t;ints_unroll_temp'];%满足基准后，将该行间隙数据存储
            end
            if ints_unroll_temp_1 < -1*E_fiber
                disp(['错误层级：', num2str(i)]);
                disp(['报错间隙：', num2str(ints_unroll_temp_1)]);
                disp(['报错匝数：', num2str(j)]);
                error_level=[error_level,i];
                %tabulate(error_level);
                break%跳出
            end
            ints_mat{i + 1,1} = ints_unroll_t;%存储间隙矩阵
            if i == 1
                coc_mat{i,1} = (R_mat{i,1})/2;
            end
            delta_d_temp_all = delta_d_temp_all(:,1:N_bot_cal);
            coc_mat{i + 1,1}=delta_d_temp_all;%存储纤芯距骨架表面距离
            
            if reconstructed_signal == 1
                disp('Data Reconstructed ?');
                turns = sum(N_bot_sat(1:(i + 1)));%发出重构信号后已用光纤长度
                samp_pos = samp_pos(1:turns*N_samp+1);
                samp_pos_temp = samp_pos_temp(1:turns+1);
                reconstructed_sg = i + 1;
                [diam_mat, samp_pos, samp_pos_temp] = reconstructed_fiber_test(i, samp_pos_control, layer, H_fiber_save, R_fiber, A_base, ints_bot, R_base, H_decr,...
                    samp_pos_temp, samp_pos, N_samp, reconstructed_sg,...
                    N_bot_sat, T_base_left, T_base_left_ex, T_base_right, conf_int, L_fiber_sat_save, L_fiber);
            end
            
        end
        if i == layer - 2
            sat = 1;%当推测完所有层级后，更换sat评判参数终止推测
        end
    end
end
