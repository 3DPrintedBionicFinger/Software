toWrite=35.6

fopen(s);
pause(3);
for i = 10:60
    fwrite(s,uint8(i));
    fscanf(s)
    pause(0.5);
end
fclose(s);