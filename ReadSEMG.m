% a = arduino('COM5','Nano');
a_pin = 'A1';
mode='AnalogInput';
fs = 500.0;   % sampling frequency (samplings per second)
mt = 1;  % time for measurements in seconds
Voltage(1:10000)=0;
TimeSampled(1:10000)=0;
SampleIndex=1;
configurePin(a,a_pin,mode)
sEMGSampleTimer=0;
tic;% start timer
while toc < mt %loop while the timer has not run mt seconds

    if((toc-sEMGSampleTimer)>(1/fs))
        Voltage(SampleIndex)=readVoltage(a,a_pin);
        sEMGSampleTimer=toc;
        TimeSampled(SampleIndex)=sEMGSampleTimer;
        SampleIndex=SampleIndex+1;
    end
    % wait for appropriate time for next measurement
    

    % beep every second
%     if (ceil(toc) > last_beep)
%         beep(); % don't know if this one exist, check docs
%         last_beep = ceil(toc);
%     end
end