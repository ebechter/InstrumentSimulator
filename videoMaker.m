figure(1)
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
for ind = 6202:-20:1
    imagesc(PSF(:,:,ind))
    drawnow
    F(ind) = getframe(gcf); 
    writeVideo(vidfile,F(ind));
end
close(vidfile)