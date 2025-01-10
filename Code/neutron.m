%一维球中子输运蒙卡模拟，计算keff
clc;clear;
%常量，变量，初始化
N = 4E3;%模拟的中子数
R = 145.5436;%裸堆半径，cm
Sigma_t = 0.05;%总截面
Sigma_s = 0.03;%散射截面
Sigma_f = 0.02;%裂变截面
nu_value = 1.125;%一次裂变产生的中子数
P_s = Sigma_s/Sigma_t;%散射发生的概率

%用一个矩阵来存储每一代中子的数据
S = zeros(7,4000);%第 1 行中子序数；第2 3 4 行位置；第5 6 7 行方向
%4000个中子初始化
S(1,:) = reshape(1:4000,1,4000);
S(7,:) = 1;
S_tobe = zeros(7,4000);%存储裂变后的 准下一代 中子
i = 0;%准下一代中子矩阵 列数,在每模拟完一代后记得重置为0
keff = zeros(200,1);%存储每一代后计算得到的keff
generation = 0;%迭代次数
nu_yiba = zeros(200,1);


while true
    %对4000个中子逐一跟踪
    for k = 1:4000
        %单个中子模拟过程
        while true
            %单个中子模拟，跳出循环的历史终止条件为发生裂变或逃出球体
            %确定下一个碰撞点位置
            rho = exponential_random_samples(1,1);%自由程指数分布抽样
            L = rho / Sigma_t;%输运距离
            S(2,k) = S(2,k) + S(5,k)*L;
            S(3,k) = S(3,k) + S(6,k)*L;
            S(4,k) = S(4,k) + S(7,k)*L;
            %判断是否超出球体边界
            if S(2,k)^2 + S(3,k)^2 + S(4,k)^2 > R
                break;
            end
            %抽取随机数，判断反应类型
            kesi_1 = rand;
            if kesi_1 > P_s
                %发生裂变反应
                kesi_2 = rand;
                if kesi_2 <= (nu_value - floor(nu_value))
                    nu_i = floor(nu_value) + 1;
                else
                    nu_i = floor(nu_value);
                end
                i = i + 1;
                S_tobe(1,i) = nu_i;%记录裂变后的中子数目
                %记录发生裂变反应时的位置
                S_tobe(2,i) = S(2,k);
                S_tobe(3,i) = S(3,k);
                S_tobe(4,i) = S(4,k);
                %结束此中子的跟踪历史
                break;
            else
                %发生散射反应
                %确定散射后方向角
                [S(5,k),S(6,k),S(7,k)] = omega_calculate(S,k);
            end
        end
    end
    generation = generation + 1;%模拟的中子代数加一
    %得到 准下一代中子 后计算keff
    nu_all = sum(S_tobe(1,:));
    nu_yiba(generation) = nu_all;
    %算出keff
    %keff(generation) = nu_all / 4000;%因为每一代的中子数量都是4000
    if generation > 1
        keff(generation) = nu_yiba(generation)/nu_yiba(generation-1);
    end
    %判断keff是否在一段时间内收敛到1，是则结束模拟，否则按权重抽样生成下一代初始中子
    %也可以只模拟200代，结束后再看这200代是否收敛
    if generation >=200
        break;
    end
    S = zeros(7,4000);
    S(1,:) = reshape(1:4000,1,4000);
    %根据 准下一代中子矩阵 抽样
    t = 0;
    for m = 1:i
        for new_num = 1:floor(4000*S_tobe(1,m)/nu_all)
            t = t + new_num;
            S(2,t) = S_tobe(2,i);
            S(3,t) = S_tobe(3,i);
            S(4,t) = S_tobe(4,i);
            %方向各向同性抽样，同时认为中子能量不变，故和裂变前速度无关
            [S(5,t),S(6,t),S(7,t)] = omega_new();
            if t == 4000
                break;
            end
        end
        if t == 4000
            break;
        end
    end
    i = 0;%重置 准下一代中子数据列数
    S_tobe = zeros(size(S_tobe));%重置 准下一代中子矩阵
end

figure(1)
plot(1:200,keff);
title('有效增殖系数变化');
xlabel('中子代数');
ylabel('有效增殖系数keff');
ylim([0,max(keff)]);
grid on;

figure(2)
plot(1:200,nu_yiba);
title('裂变后中子数量变化');
xlabel('中子代数');
ylabel('裂变后产生的中子数量');
ylim([0,max(nu_yiba)]);
disp(mean(keff(20:200)));
grid on;

