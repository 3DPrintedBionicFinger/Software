%s=serial('COM5','BaudRate',128000);


mt=10;

SerialData={};
SerialIndex=1;
fopen(s);
SerialData=[SerialData
        fgetl(s)];
pause(3);
tic;% start timer
while toc < mt %loop while the timer has not run mt seconds
    SerialData=[SerialData
        fgetl(s)];
        
end
fclose(s);
length(SerialData)/mt