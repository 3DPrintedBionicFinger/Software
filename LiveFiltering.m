%s=serial('COM5','BaudRate',128000);

SRate=500;

f2=SRate*(0:(FFTsize-1))/FFTsize;
loopCount=1;

startTime=0;

freqMax=160;
freqMin=60;
FFTSize=128;
toWrite=0;
mt=100;

d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',29,'HalfPowerFrequency2',63, ...
               'DesignMethod','butter','SampleRate',SRate);
SerialData={};
freqIndex=1;
fopen(s);
pause(3);%needs to set up the port
fgetl(s);%should getrid of any bad data at the start
fgetl(s);
fgetl(s);
fgetl(s);
fgetl(s);
fgetl(s);
fgetl(s);
fgetl(s);
fgetl(s);
tic;% start timer
while toc < mt %loop while the timer has not run mt seconds
    SerialData=[SerialData
        fgetl(s)];
%     if (str2double(SerialData(end))>680) || (str2double(SerialData(end))<640)
%         SerialData(end)
%     end
    if (length(SerialData)>(FFTSize/2+1))
        sampleTime=toc-startTime;
        dataToFilter=str2double(SerialData(1:FFTSize/2));%the start points help avoid band transfers at the start
        
        %filtering
        %dataToFilter=filtfilt(d,dataToFilter);%applies notch filter
        dataToFilter=dataToFilter-mean(dataToFilter);%applies DC filter
        wave=modwt(dataToFilter','db2',2);
        
        %FFT
        waveFFT=abs(fft(wave(2,:),FFTSize));
        
        
        [mag,index]=max(waveFFT);
        freq(freqIndex)=(index-1)*SRate/(FFTSize);
        freqIndex=mod(freqIndex,4)+1;
        
        avgFreq=sum(freq(1:4))/4
        
        
        
        %stuff that is not strictly nessicariry
        dataFFT=abs(fft(dataToFilter,FFTSize));
        %figure;
        subplot(2,1,1);
        plot(f2(1:end/2),waveFFT(1:(FFTSize/2)));
        hold on
        plot(f2(1:end/2),dataFFT(1:(FFTSize/2)));
        hold off
        subplot(2,1,2);
        plot(wave(2,:))
        hold on
        plot(dataToFilter)
        hold off
        loopCount=loopCount+1;
            pause(0.05);
        
        
        if (avgFreq>freqMax)
            avgFreq=freqMax;
        end
        if (avgFreq<freqMin)
            avgFreq=freqMin;
        end
        toWrite=uint8(((avgFreq-freqMin)/(freqMax-freqMin))*255);
        
        fwrite(s,toWrite);
        
        SerialData=SerialData(FFTsize/2:end);
        startTime=toc;
    end 
end
fclose(s);
loopCount