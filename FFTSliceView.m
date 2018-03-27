FFTsize=128;

%data=csvread("MikeFlex-10S-675Hz");
%data=csvread("MikeRest-20S-731Hz");
data=sEMG10;

% sampleTimes(1:999)=0;
% for x= 1:1000
%     sampleTimes(x)=(data(x+1,2)-data(x,2));
% end
% temp=mean(sampleTimes);
% temp=temp/(1000000);
% SRate=round(1/temp);

SRate=661;
f2=SRate*(0:(FFTsize-1))/FFTsize;
clear freq
freq(1:floor(length(data)/FFTsize))=0;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',29,'HalfPowerFrequency2',67, ...
               'DesignMethod','butter','SampleRate',SRate);
%data=filtfilt(d,data);
for m = 1:floor(length(data)/FFTsize)
    pause(1)
%     figure
    tic;
    data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1)=filtfilt(d,data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1));
    DC=mean(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1));
    for n = (FFTsize+1)*(m-1)+1:(FFTsize+1)*(m)
        data(n,1)=data(n,1)-DC;
    end
    wave=modwt(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1)','db2',2);
%     DC=mean(data(FFTsize*m:FFTsize*(m+1),1));
%     data(FFTsize*m:FFTsize*(m+1))=(FFTsize*m:FFTsize*(m+1))-DC;

    Y5=fft(wave(2,:),FFTsize);
    toc
    Y6=fft(data((FFTsize+1)*(m-1)+1:(FFTsize+1)*(m),1),FFTsize);
   % hold on
    %plotData=figure;
    subplot(2,1,1);
    plot(f2(1:end/2),abs(Y5(1:(FFTsize/2))));
    [MAG,IND]=max(abs(Y5(1:(FFTsize/2))));
    freq(m)=(IND-1)*SRate/(FFTsize);
    (IND-1)*SRate/(FFTsize)
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
rollAvg=movmean(freq,4)'
