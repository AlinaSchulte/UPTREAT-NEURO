function stim = ud_chirps(dur,fs,fl,fu,len)

 t = 0:1/fs:len-1/fs;

 fmod = .5+.5*sin((2/dur)*pi*t+3*pi/2); 
 fmod = fmod*(fu-fl)+fl; 
 phi = 2 * pi * cumsum(fmod) / fs; 
 stim = sin(phi);
 end