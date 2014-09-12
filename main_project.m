av=load('testsample1.mfcc');                                 %  loading the acoustic vector into a
global speaker_pos
len_av=length(av);                                      %  length of total accoustic vectors
if (rem(len_av,10)>4)
    no_of_win=((len_av-rem(len_av,5)-25)/10);          %  number of windows
else
    no_of_win=((len_av-rem(len_av,5)-30)/10);
end

no_of_win_2=((len_av-rem(len_av,10)-50)/10);                                %  number of windows of double size
win=zeros(no_of_win,25,13);                              %  array of windows 
win_2=zeros(no_of_win_2,50,13);                         %  array of windows of double the size


%  to create a array of windows which are of 5x13...........................
dumv=zeros(25,13);                                       %  dummy vector which will help in transferring data into win
for i=1:no_of_win
    for j=1:25
        dumv(j,:)=av((i-1).*10+j,:);
    end
    win(i,:,:)=dumv(:,:); 
end

%  to create a array of windows of double the size which are of 5x13........
dumv_2=zeros(50,13);                                    %  dummy vector which will help in transferring data into win_2
for i=1:no_of_win_2
    for j=1:50
        dumv_2(j,:)=av((i-1).*10+j,:);
    end
    win_2(i,:,:)=dumv_2(:,:); 
end


mean_win=zeros(no_of_win,13);                           %  array of means of each window
mean_win_2=zeros(no_of_win_2,13);                       %  array of means of each window of double the size
cov_win=zeros(no_of_win,13,13);                         %  array of covariance of each window
cov_win_2=zeros(no_of_win_2,13,13);                     %  array of covariance of each window of double the size

%  to calculate array of means and covariance of each window ...............
dumc=zeros(25,13);
for i=1:no_of_win
    mean_win(i,:)=mean(win(i,:,:));
    dumc(:,:)=win(i,:,:);                               %  dummy to calculate covariance
    cov_win(i,:,:)=cov(dumc(:,:));
end

%  to calculate array of means and covariance of each window which are of double the size 
dumc_2=zeros(50,13);
for i=1:no_of_win_2
    mean_win_2(i,:)=mean(win_2(i,:,:));
    dumc_2(:,:)=win_2(i,:,:);                           %  dummy to calculate covariance
    cov_win_2(i,:,:)=cov(dumc_2(:,:));                  
end

win_pdf=zeros(no_of_win,25);
win_pdf_2=zeros(no_of_win_2,50);
%  creates the p-d-f with values at respective numbers given ..............

dum_1_1=zeros(25,13);
dum_2_1=zeros(1,13);
dum_3_1=zeros(13,13);

for i =1:no_of_win
    dum_1_1(:,:)=win(i,:,:);
    dum_2_1(:,:)=mean_win(i,:);
    dum_3_1(:,:)=cov_win(i,:,:);
    win_pdf(i,:)=mvnpdf(dum_1_1,dum_2_1,dum_3_1);
end

%  creates the p-d-f for windows of double the size with values at respective numbers given 

dum_1_2=zeros(50,13);
dum_2_2=zeros(1,13);
dum_3_2=zeros(13,13);

for i =1:no_of_win_2
    
    dum_1_2(:,:)=win_2(i,:,:);
    dum_2_2(:,:)=mean_win_2(i,:);
    dum_3_2(:,:)=cov_win_2(i,:,:);
    win_pdf_2(i,:)=mvnpdf(dum_1_2,dum_2_2,dum_3_2);
end

win_lik=ones(1,no_of_win);                             %  array of how likely the window fits into the distribution                        
win_lik_2=ones(1,no_of_win_2);                           %  array of how likely the window of double the size fits into the distribution                        

%  calculates how likely the window fits into the distribution.............
for i=1:no_of_win
    for j=1:5
        win_lik(i)=win_lik(i).*win_pdf(i,j);        
    end
end

%  calculates how likely the window fits into the distribution.............
for i=1:no_of_win_2
    for j=1:10
        win_lik_2(i)=win_lik_2(i).*win_pdf_2(i,j);  
    end
end

%  to calculate the likelihood ratio for each window........................
lambda_win=zeros(1,no_of_win_2);

for i=1:no_of_win_2
    lambda_win(i)=(win_lik_2(i)/(win_lik(i).*win_lik(i+1)));                                    
end


%  to calculate the distance for each window...............................
dist_win=zeros(1,no_of_win_2);

for i=1:no_of_win_2
    dist_win(i)=-log(lambda_win(i));
end

%plot(dist_win);


%  to filter the dist_win by considering mean of 10 windows

filt_poi=zeros(1,(length(dist_win)-rem(length(dist_win),10)) /10);

for i=1:length(filt_poi)
    for j=1:10
        filt_poi(i)=(dist_win((i-1)*10+j)+filt_poi(i));
    end
    filt_poi(i)=filt_poi(i)/10;
end

%plot(filt_poi);

len_filt=length(filt_poi);
ch_pos_dis_h=zeros(1,len_filt);
k=1;
for i=1:len_filt
    if(i==1)
        if(filt_poi(i+1)<filt_poi(i))
            ch_pos_dis_h(k)=i;
            k=k+1;
        end
    elseif(i==len_filt)
        if(filt_poi(i-1)<filt_poi(i))
            ch_pos_dis_h(k)=i;
             k=k+1;
        end
    elseif(filt_poi(i)>filt_poi(i-1) && filt_poi(i)>filt_poi(i+1))
            ch_pos_dis_h(k)=i;
            k=k+1;
    end
end

ch_pos_dis=zeros(1,k-1);
for i=1:k-1
    ch_pos_dis(i)=((10.*ch_pos_dis_h(i))-5).*10+25;
end



% second pass into bic for futher refining.......................

lambda=3.5;

s1=1;
dc=1;
ch_pos_bic_h=zeros(1,k-1);
e1=ch_pos_dis(1);
s2=e1+1;
e2=ch_pos_dis(2);
for d=1:12
    ch_pos_bic_h(dc)=s2;
    del_val=bic_func(s1,e1,s2,e2,av,lambda);
    if(del_val<0)
        
        %a=1
        s1=s2;
        e1=e2;
        s2=e1+1;
        if(d==k-1)
            e2=length(av)-10;
        else
        
            if(d==k)
                break
            end
            e2=ch_pos_dis(d+1);
        end
        dc=dc+1;
        
    else
        e1=e2;
        s2=e1+1;
       
        if(d==k-1)
            e2=length(av)-10;
        else
            if(d==k)
                break
            end
            e2=ch_pos_dis(d+1);
        end
    end
    
end

ch_pos_bic=zeros(1,dc);
for i=1:dc
    ch_pos_bic(i)=ch_pos_bic_h(i);
end

%plot(ch_pos_bic);
%..........clustering the speech change points

%..........assumed to have dc+1 different speakers so each one of them are
%..........assigned a starting and ending points

sp_pos=zeros(dc+1,dc+1,2);                            %position of start and end points of speaker in window id's
sp_pos(1,1,1)=1;

for i=1:dc
    sp_pos(i,1,2)=ch_pos_bic(i)-1;
    sp_pos(i+1,1,1)=ch_pos_bic(i);
end
sp_pos(dc+1,1,2)=length(av(:,:));

clu_h=zeros(dc+1,dc+1);

for i=1:dc
    for j=i+2:dc
        if(i==1)
            s1=1;
            e1=ch_pos_bic(i);
            s2=ch_pos_bic(j-1)+1;
            e2=ch_pos_bic(j);
        else
            s1=ch_pos_bic(i-1)+1;
            e1=ch_pos_bic(i);
            s2=ch_pos_bic(j-1)+1;
            e2=ch_pos_bic(j);
        end
        clu_h(i,j)=bic_func(s1,e1,s2,e2,av,lambda);    
    end
end

%sp_pos=zeros(dc+1,dc+1,2);
for i=1:dc+1
    k=2;
    for j=i+2:dc
        if(clu_h(i,j)>0)
            if(sp_pos(i,1,1)~=0)
                sp_pos(i,k,1)=ch_pos_bic(j-1)+1;
                sp_pos(i,k,2)=ch_pos_bic(j);
            end
            sp_pos(j,1,1)=0;
            sp_pos(j,1,2)=0;
            k=k+1;
        end
    end
end

speaker_position_h=zeros(dc+1,dc+1,2);
for i=1:dc+1
    for j=1:dc+1
        if(sp_pos(i,j,1)~=0)
         %   speaker_positions(i,j,1)=1;
            speaker_position_h(i,j,1)=sp_pos(i,j,1).*480;
        end
        if(sp_pos(i,j,2)~=0)
         %   speaker_positions(i,j,1)=1;
            speaker_position_h(i,j,2)=sp_pos(i,j,2).*480;
        end
    end
end
 
N=0;
for i=1:dc+1
    if(speaker_position_h(i,1,1)~=0)
        N=N+1;
    end
end

speaker_pos=zeros(N,dc+1,2);
sp=zeros(1,N);
k=1;
for i=1:dc+1
    if(speaker_position_h(i,1,1)~=0)
        j=1;
        while(speaker_position_h(i,j,1)~=0)
            speaker_pos(k,j,1)=speaker_position_h(i,j,1);
            speaker_pos(k,j,2)=speaker_position_h(i,j,2);
            sp(k)=sp(k)+(speaker_position_h(i,j,2)-speaker_position_h(i,j,1)+1);
            j=j+1;
        end
        k=k+1;
    end
end

max=sp(1);
for i=2:N
    if(sp(i)>max)
        max=sp(i);
    end
end


%{
 samp_data=wavread('sample4.wav');
 player=audioplayer(samp_data,48000);
 %Speaker1
 
 playblocking(player,[speaker_pos(1,1,1) speaker_pos(1,1,2)]); 
 playblocking(player,[speaker_pos(1,2,1) speaker_pos(1,2,2)]); 
 playblocking(player,[speaker_pos(1,3,1) speaker_pos(1,3,2)]);
 play(player,[speaker_pos(1,4,1) speaker_pos(1,4,2)]);
 
 %pause(player);
 %Speaker2
 resume(player);
 playblocking(player,[speaker_pos(2,1,1) speaker_pos(2,1,2)]); 
 play(player,[speaker_pos(2,2,1) speaker_pos(2,2,2)]); 
 %playblocking(player,[speaker_pos(2,3,1) speaker_pos(2,3,2)]);
 %play(player,[speaker_pos(1,4,1) speaker_pos(1,4,2)]);
  %}
%stgui();

 
                    