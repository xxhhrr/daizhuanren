%�û���Ŵ��㷨�Ż���������������
%�ο����ף�MIMO�״�����������Ƽ��źŴ����о�
%�����ƽ����ֵ�԰�Ϊ-16.7dB���ó�����ɴﵽ-13.8dB���ң�ƽ������ط�ֵ��ƽ������ֵ��ࡣ
%���ӿƼ���ѧ2018��׿Խ�����Ȼ
close all;
N = 40;      %�źų���
L = 4;       %�źŸ���
M = 4;       %M����
mutationProb=0.001;   %�������
Num = 200;   %��Ⱥ����
Q = zeros(Num,4);
w = ones(1,4);
d = zeros(Num,4);
E = zeros(1,Num);
E0 = 1;    %����Ŀ�����Ӧ��ֵ
%������һ����Ⱥ
st_pop = zeros(Num*L,N+1);
st_pop = pi/2*round(rand(Num*L,N)*4-0.50000000001);
times = 5000;   %���ѡ�����
t = 0;
cost_value = zeros(1,Num*L-3);

for ii = 1:Num
    Q(ii,:) = cost(st_pop(4*ii-3:4*ii,:),L,N,w);
    E(1,ii) = Q(ii,1)*w(1) + Q(ii,2)*w(2) + Q(ii,3)*w(3) + Q(ii,4)*w(4);
end
A = min(E);
[~,n] = find(A==min(E));        %�ҳ���õĸ���
w(2) = Q(n,1)/Q(n,2);           %����w
w(3) = (Q(n,1)/Q(n,3))^2;
w(4) = (Q(n,1)/Q(n,4))^2;
for ii = 1:Num
    E(1,ii) = Q(ii,1)*w(1) + Q(ii,2)*w(2) + Q(ii,3)*w(3) + Q(ii,4)*w(4);    
end
%%���̶�ѡ���㷨
flag = 1;
st = 1;
pop2 = zeros(1,N);
P = zeros(1,N);
cha = 0;
for t = 1:times
    if min(E) > E0 || cha > 0.1
        Et = min(E);
        E1 = max(E) - E;    %���̶�ѡ����
        Size = size(E1);
        for ii = 1:Size(2)
            P(ii) = E1(ii)/sum(E1)*100;  %�Ŵ�����
            p = rand(1,1);
            if P(ii) > p
                %�µ���Ⱥ
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
        %����
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
            %����������Ӧ�Ⱥ���
            Q2 = zeros(Length(1),4);
            E = zeros(1,Length(1)/4);
            for ii = 1:4:Length(1)
                Q2(ii,:) = cost(pop2(ii:ii+3,:),L,N,w);     
                E(1,(ii+3)/4) = Q2(ii,1) * w(1) +Q2(ii,2) * w(2) + Q2(ii,3) * w(3) + Q2(ii,4) * w(4); %�ڶ����Ĵ��ۺ�
            end
            st_pop = pop2;
            pop2 = zeros(1,N);
            st = 1;
            flag = 1;
        end
            cha(t) = norm(Et-min(E));
    end
end
%�ڶ����Ż�-����������
A = min(E);
[m,n] = find(A==min(E));
tar_pop = st_pop(4*n-3:4*n,:);   %��һ�������Ż������   
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
subplot(221),plot(20*log10(Q8(1,:)));title('����1������غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(222),plot(20*log10(Q8(2,:)));title('����2������غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(223),plot(20*log10(Q8(3,:)));title('����3������غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(224),plot(20*log10(Q8(4,:)));title('����4������غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
figure(2)
subplot(321),plot(20*log10(Q9(1,:)));title('����1������2�Ļ���غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(322),plot(20*log10(Q9(2,:)));title('����1������3�Ļ���غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(323),plot(20*log10(Q9(3,:)));title('����1������4�Ļ���غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(324),plot(20*log10(Q9(4,:)));title('����2������3�Ļ���غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(325),plot(20*log10(Q9(5,:)));title('����2������4�Ļ���غ���');ylabel('����(dB)');
axis([-inf,inf,-30,0]);grid on;
subplot(326),plot(20*log10(Q9(6,:)));title('����3������4�Ļ���غ���');ylabel('����(dB)');
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