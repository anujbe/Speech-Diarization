function speaker1() 
global speaker_pos samp_data
player=audioplayer(samp_data,48000);
%samp_data = 0.99*samp_data/max(abs(samp_data));
       
        playblocking(player,[speaker_pos(1,1,1) speaker_pos(1,1,2)]);
        playblocking(player,[speaker_pos(1,2,1) speaker_pos(1,2,2)]);
        playblocking(player,[speaker_pos(1,3,1) speaker_pos(1,3,2)]);
        %play(player,[speaker_pos(1,4,1) speaker_pos(1,4,2)]);
        %playblocking(player,[speaker_pos(1,5,1) speaker_pos(1,5,2)]);
 
 