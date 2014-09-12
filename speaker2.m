function speaker2() 
global speaker_pos samp_data
 
 player=audioplayer(samp_data,48000);
 
 playblocking(player,[speaker_pos(2,1,1) speaker_pos(2,1,2)]); 
 playblocking(player,[speaker_pos(2,2,1) speaker_pos(2,2,2)])
 %play(player,[speaker_pos(2,3,1) speaker_pos(2,3,2)]); 
  