function Q = cost(pop1,L,N,w)
B = 0;
Q1 = 0;
Q2 = 0;
pop = pop1;
for l = 1:L
    for k = -N+1:N-1
        if k == 0
            continue
        end
          B = max(B,abs(Cal_Cor(pop(l,:),pop(1,:),k,N)));
    end
    Q1 = Q1 + B;
end
B = 0;
x = 0;
y = 0;
for p=1:L-1
    for q=p+1:L
        for k=(1-N):(N-1)
            if(p==q)
                continue;
            end
            sp = pop(p,:);
            sq = pop(q,:);
            x = Cal_Cor(sp,sq,k,N);
            B = max(norm(x),B);
        end
        y = y + B;
    end
    Q2 = Q2 + y;
end
x = 0;
y = 0;
B = 0;
Q3 = 0;
for l = 1:L
    for k = 1:N-1
        x = Cal_Cor(pop(l,:),pop(1,:),k,N);
        B = (norm(x))^2 + B;
    end
    Q3 = Q3 + B;
end
B = 0;
x = 0;
y = 0;
Q4 = 0;
for p=1:L-1
    for q=p+1:L
        for k=(1-N):(N-1)
            sp = pop(p,:);
            sq = pop(q,:);
            x = Cal_Cor(sp,sq,k,N);
            B = (norm(x))^2 + B;
        end
        y = y + B;
    end
    Q4 = Q4 + y;
end
E = w(1) * Q1 + w(2) * Q2 + w(3) * Q3 + w(4) * Q4;
Q = [Q1 Q2 Q3 Q4];
