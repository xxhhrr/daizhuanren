%用混合遗传算法优化产生正交多项码
%参考文献：MIMO雷达正交波形设计及信号处理研究
%理想的平均峰值旁瓣为-16.7dB，该程序仅可达到-13.8dB左右，平均互相关峰值电平与理论值差不多。
%电子科技大学2018级卓越班徐浩然
close all;
N = 40;      %信号长度
L = 4;       %信号个数
M = 4;       %M相码
mutationProb=0.001;   %变异概率
Num = 200;   %种群数量
Q = zeros(Num,4);
w = ones(1,4);
d = zeros(Num,4);
E = zeros(1,Num);
E0 = 1;    %满足目标的适应度值
%产生第一代种群
st_pop = zeros(Num*L,N+1);
st_pop = pi/2*round(rand(Num*L,N)*4-0.50000000001);
times = 5000;   %最大选择次数
t = 0;
cost_value = zeros(1,Num*L-3);

for ii = 1:Num
    Q(ii,:) = cost(st_pop(4*ii-3:4*ii,:),L,N,w);
    E(1,ii) = Q(ii,1)*w(1) + Q(ii,2)*w(2) + Q(ii,3)*w(3) + Q(ii,4)*w(4);
end
A = min(E);
[~,n] = find(A==min(E));        %找出最好的个体
w(2) = Q(n,1)/Q(n,2);           %调整w
w(3) = (Q(n,1)/Q(n,3))^2;
w(4) = (Q(n,1)/Q(n,4))^2;
for ii = 1:Num
    E(1,ii) = Q(ii,1)*w(1) + Q(ii,2)*w(2) + Q(ii,3)*w(3) + Q(ii,4)*w(4);    
end
%%轮盘赌选择算法
flag = 1;
st = 1;
pop2 = zeros(1,N);
P = zeros(1,N);
cha = 0;
for t = 1:times
    if min(E) > E0 || cha > 0.1
        Et = min(E);
        E1 = max(E) - E;    %轮盘赌选择函数
        Size = size(E1);
        for ii = 1:Size(2)
            P(ii) = E1(ii)/sum(E1)*100;  %遗传概率
            p = rand(1,1);
            if P(ii) > p
                %新的种群
                x = st_pop(4*ii-3,5:35);
                st_pop(4*ii-3,5:35) = st_pop((4*ii-2),5:35);
                st_pop(4*ii-2,5:35) = x;
                x = st_pop(4*ii-1,5:35);
                st_pop(4*ii-1,5:35) = st_pop((4*ii),5:35);
                st_pop(4*ii,5:35) = x;
                pop2 = [pop2 
                        st_pop(4*ii-3:4*ii,:)];
                if st == 1
                    pop2 = pop2(2:5,:);
                    st = 0;
                end
                flag = 0;
            end
        end
        %变异
        if flag == 0
            Length = size(pop2(:,1));
            for l = 1:Length(1)
                for n = 1:N
                    Prob = rand(1,1);
                    if Prob < mutationProb
                        pop2(l,n) = pi*round(rand(1,1)*3-0.5);
                    end
                end
            end
            %重新评估适应度函数
            Q2 = zeros(Length(1),4);
            E = zeros(1,Length(1)/4);
            for ii = 1:4:Length(1)
                Q2(ii,:) = cost(pop2(ii:ii+3,:),L,N,w);     
                E(1,(ii+3)/4) = Q2(ii,1) * w(1) +Q2(ii,2) * w(2) + Q2(ii,3) * w(3) + Q2(ii,4) * w(4); %第二代的代价函
            end
            st_pop = pop2;
            pop2 = zeros(1,N);
            st = 1;
            flag = 1;
        end
            cha(t) = norm(Et-min(E));
    end
end
%第二步优化-叠代码搜索
A = min(E);
[m,n] = find(A==min(E));
tar_pop = st_pop(4*n-3:4*n,:);   %第一个部分优化的输出   
for ii = 1:L
    for n = 1:N
        [A,tar_pop(ii,n)] = diedai(tar_pop,ii,n,L,N,w);
    end
end
B = 0;
Q8 = zeros(L,2*N-1);
for l = 1:L
    for k = -N+1:N-1
        B(k+N) = abs(Cal_Cor(tar_pop(l,:),tar_pop(l,:),k,N));
    end
    Q8(l,:) = B;
end
B = 0;
Q9 = zeros(6,2*N-1);
y = 1;
for p=1:L-1
    for q=p+1:L
        for k=(1-N):(N-1)
            if(p==q)
                continue;
            end
            sp = tar_pop(p,:);
            sq = tar_pop(q,:);
            x = Cal_Cor(sp,sq,k,N);
            B(k+N) = abs(x);
        end
        if(p~=q)
            Q9(y,:) = B(:);
            y = y + 1;
        end
    end
end    
figure(1)
subplot(221),plot(20*log10(Q8(1,:)));title('序列1的自相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(222),plot(20*log10(Q8(2,:)));title('序列2的自相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(223),plot(20*log10(Q8(3,:)));title('序列3的自相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(224),plot(20*log10(Q8(4,:)));title('序列4的自相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
figure(2)
subplot(321),plot(20*log10(Q9(1,:)));title('序列1和序列2的互相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(322),plot(20*log10(Q9(2,:)));title('序列1和序列3的互相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(323),plot(20*log10(Q9(3,:)));title('序列1和序列4的互相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(324),plot(20*log10(Q9(4,:)));title('序列2和序列3的互相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(325),plot(20*log10(Q9(5,:)));title('序列2和序列4的互相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(326),plot(20*log10(Q9(6,:)));title('序列3和序列4的互相关函数');ylabel('幅度(dB)');
axis([-inf,inf,-30,0]);grid on;
aa = zeros(4,4);
aa(1,1) = max(Q8(1,:));
aa(1,2) = max(Q9(1,:));
aa(1,3) = max(Q9(2,:));
aa(1,4) = max(Q9(3,:));

aa(2,1) = max(Q9(1,:));
aa(2,2) = max(Q8(2,:));
aa(2,3) = max(Q9(4,:));
aa(2,4) = max(Q9(5,:));

aa(3,1) = max(Q9(2,:));
aa(3,2) = max(Q9(4,:));
aa(3,3) = max(Q8(3,:));
aa(3,4) = max(Q9(6,:));

aa(4,1) = max(Q9(3,:));
aa(4,2) = max(Q9(5,:));
aa(4,3) = max(Q9(6,:));
aa(4,4) = max(Q8(4,:));