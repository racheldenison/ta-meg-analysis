% rd_simIntermodulation.m

dt = .001;
f1 = 30;
f2 = 40;
t = (1:1000)*dt;
f = 0:999;
 
s1 = sin(t*2*pi*f1);
s2 = sin(t*2*pi*f2);
% s3 = s1.*s2;
% s3 = (s1+s2)/2+s1.*s2;
% s3 = (s1+s2)/2+s1.*(2*s2);
% s3 = (s1+s2)/2+s1./s2;
% s3 = (s1+s2)/2+s1.^s2;
% s3 = (s1+s2)/2+s2.^s1;
s3 = s1.^2 + s2.^2;

 
figure; 
subplot(2,1,1)
plot(t,s1, t, s2, t, s3);
legend('S1', 'S2', 'S3')
 
subplot(2,1,2)
plot(f,abs(fft(s1)), f,abs(fft(s2)), f,abs(fft(s3)))
xlim([0 80])
% legend('S1', 'S2', 'S1.*S2')
% legend('S1', 'S2', 'S1+S2+S1.*S2')
% legend('S1', 'S2', 'S1+S2+S1.*(2*S2)')
% legend('S1', 'S2', 'S1+S2+S1./S2')
% legend('S1', 'S2', 'S1+S2+S1.^S2')
legend('S1', 'S2', 'S1+S2+S2.^S1')


% % demo of convolution in the frequency domain
% f = -200:200;
% f1 = zeros(size(f));
% f1(f==30) = 1;
% f1(f==-30) = 1;
% f2 = zeros(size(f));
% f2(f==40) = 1;
% f2(f==-40) = 1;
% fc = conv(f1, f2);
% 
% figure
% hold on
% plot(f, f1, f, f2)
% plot(f(1)*2:f(end)*2, fc,'r')
% xlim([f(1) f(end)])
% legend('f1','f2','conv(f1,f2)')

