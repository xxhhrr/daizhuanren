function A = Cal_Cor(sp,sq,k,N)
% Calculate the cross correlation function
% sp--(in)the phase code of one signal,1*N
% sq--(in)the phase code of another signal,1*N
% k--(in)delay,(-N+1)<=k<=(N-1)
% N--(in)length of the signal
% A--(out)the cross correlation function?
% attention: if(p==q) then it is the autocorrelative function.
if((k<N)&&(k>=0))
    A = sum(exp(1j*(sp(1,1:N-k)-sq(1,k+1:N))))/N;
elseif((k<0)&&(k>-N))
    A = sum(exp(1j*(sp(1,1-k:N)-sq(1,1:N+k))))/N;
else
    disp('Error!');
end