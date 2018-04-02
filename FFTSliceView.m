FFTsize=128;

%data=csvread("MikeFlex-10S-675Hz");
%data=csvread("MikeRest-20S-731Hz");
data=downsample(double(emg44(:,1)),4);
%data=sEMG11;

% sampleTimes(1:999)=0;
% for x= 1:1000
%     sampleTimes(x)=(data(x+1,2)-data(x,2));
% end
% temp=mean(sampleTimes);
% temp=temp/(1000000);
% SRate=round(1/temp);

SRate=500;
f2=SRate*(0:(FFTsize-1))/FFTsize;
clear freq
freqToGraph=0;
freq(1:floor(length(data)/FFTsize))=0;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',29,'HalfPowerFrequency2',67, ...
               'DesignMethod','butter','SampleRate',SRate);
%data=filtfilt(d,data);
for m = 1:floor(length(data)/FFTsize)-1
%     pause(1)
%     figure
    data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1)=filtfilt(d,data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1));
    DC=mean(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1));
    for n = (FFTsize+1)*(m-1)+1:(FFTsize+1)*(m)
        data(n,1)=data(n,1)-DC;
    end
    wave=modwt(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1)','db2',2);
%     DC=mean(data(FFTsize*m:FFTsize*(m+1),1));
%     data(FFTsize*m:FFTsize*(m+1))=(FFTsize*m:FFTsize*(m+1))-DC;

    Y5=fft(wave(2,:),FFTsize);
    Y6=fft(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1),FFTsize);
   % hold on
    %plotData=figure;
    subplot(2,1,1);
    plot(f2(1:end/2),abs(Y5(1:(FFTsize/2))));
    [MAG,IND]=max(abs(Y5(1:(FFTsize/2))));
    freq(m)=(IND-1)*SRate/(FFTsize);
    (IND-1)*SRate/(FFTsize);
    hold on
    plot(f2(1:end/2),abs(Y6(1:(FFTsize/2))));
    hold off
    subplot(2,1,2);
    plot(wave(2,:))

    hold on
    plot(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m)))
    hold off
    %saveas(plotData,sprintf('sEMG7-FFT-%d.jpeg',m));
    
end
rollAvg=movmean(freq,8);
for k=1:floor(length(data)/FFTsize)-20
    for i= 1:FFTsize
        freqToGraph(((k-1)*FFTsize)+i)=(rollAvg(k)-50)/100000;
    end
end
figure

plot(data(1:end-FFTsize*20))
hold on
plot(freqToGraph')
hold off
h=0;
h2=0;
h3=0;
for k= 1:length(freq)-1
    h(k)=freq(k+1)-freq(k);
    h2(k)=rollAvg(k+1)-rollAvg(k);
    h3(k)=rollAvg(k+1)-freq(k);
    if (h2(k)>0)
        1
    else
        0
    end
end
% figure
% plot(h);
figure
plot(h2);
% figure
% plot(h3);
