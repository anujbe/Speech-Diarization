%{
function speaker4()
global speaker_pos samp_data

 %player=audioplayer(samp_data,48000);
 
 playblocking(player,[speaker_pos(1,4,1) speaker_pos(1,4,2)]); 
 
 playblocking(player,[speaker_pos(1,2,1) speaker_pos(1,2,2)]); 
 playblocking(player,[speaker_pos(1,3,1) speaker_pos(1,3,2)]);
 play(player,[speaker_pos(1,4,1) speaker_pos(1,4,2)]);
 %}