implay('DSCF6004.MOV');
v = VideoReader('DSCF6004.MOV');
k=1;
video = zeros(450,700,3,2550,'uint8');
while hasFrame(v)
    frame = readFrame(v);
    video(:,:,:,k)= frame(401:850,651:1350,:);
    k=k+1;
    if ~mod(k,100)
        k
    end
end

aviFile = VideoWriter('newfile.avi');
open(aviFile);
writeVideo(aviFile,video);
close(aviFile);
