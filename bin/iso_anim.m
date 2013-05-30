function iso_anim(start,finish)
for i=start:finish
    gather(i)
    var_iso
    fOUT=sprintf('./data/iso%03d.jpeg',i)
    print('-djpeg',fOUT)
    close all
end