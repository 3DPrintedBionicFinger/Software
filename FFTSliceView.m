FFTsize=128;

data=sEMG;


% sampleTimes(1:999)=0;
% for x= 1:1000
%     sampleTimes(x)=(data(x+1,2)-data(x,2));
% end
% temp=mean(sampleTimes);
% temp=temp/(1000000);
% SRate=round(1/temp);

SRate=500;
f2=SRate*(0:(FFTsize-1))/FFTsize;


for m = 1:20
    pause(1)
    DC=mean(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1));
    for n = (FFTsize+1)*(m-1)+1:(FFTsize+1)*(m)
        data(n,1)=data(n,1)-DC;
    end
%     DC=mean(data(FFTsize*m:FFTsize*(m+1),1));
%     data(FFTsize*m:FFTsize*(m+1))=(FFTsize*m:FFTsize*(m+1))-DC;
    Y5=fft(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1),FFTsize);
    hold on
    %plotData=figure;
    plot(f2(1:end/2-30),abs(Y5(1:(FFTsize/2-30))));
    %saveas(plotData,sprintf('sEMG7-FFT-%d.jpeg',m));
end