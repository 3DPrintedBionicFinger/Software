FFTsize=128;
SRate=500;
f2=SRate*(0:(FFTsize-1))/FFTsize;
data=sEMG2;
data2=sEMG2(:,1);
for m = 1:1
    pause(1)
    DC=mean(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1));
    for n = (FFTsize+1)*(m-1)+1:(FFTsize+1)*(m)
        data(n,1)=data(n,1)-DC;
    end
    %wave=modwt(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1)','db45',4);
    %Y5=fft(wave(2,:),FFTsize);
    Y6=fft(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1),FFTsize);
   % hold on
    %plotData=figure;
    %plot(f2(1:end/2),abs(Y5(1:(FFTsize/2))));
    hold on
    plot(f2(1:end/2),abs(Y6(1:(FFTsize/2))));
    hold off
    %saveas(plotData,sprintf('sEMG7-FFT-%d.jpeg',m));
end
[c,l]=wavedec(data2,3,'db2');
approx = appcoef(c,l,'db2');
plot(approx)
hold on
plot(data2(:,1))
plot(abs(fft(data2)))