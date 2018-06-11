
%%%%%% IITH and PiiTech collaboration work %%%%%%%%%%
%%%%%% Below code is to find boundaries of each & every beat(1Khz),
%%%%%% followed by Feature-Extraction of all beats%%%%%%%%%%%%%%%%%%%

% function Updated_Boundary_Detection_to_PiiTech(ECG,Fs,name)

clear all
close all
clc

Fs = 122;

% load('E:\ECG_Waveforms_Matlab\PiiTech_Data\Niranjan_database_copy\mitdb_217_181_186.mat');
% ECG= signal(:,2);
% ECG=ECG(1:50000);

ECG = xlsread('mitdb_109_15_20.xls',1,'A1:A610');
%   ECG = xlsread('E:\ECG_Waveforms_Matlab\PiiTech_Data\mitdb_108_9_14.csv',1,'A1:A610');
% E:\ECG_Waveforms_Matlab\PiiTech_Data
% mitdb_230_44_49
%  ECG = xlsread('E:\ECG_Waveforms_Matlab\PiiTech_Data\Niranjan_database\Passed\mitdb_217_181_186.csv',1,'A1:A610');

% val=val';
% ECG=val;q

ECG

figure(1);
plot(ECG)
title('input');

xnew1=filtering(ECG);
xnew2=xnew1(:,1);
xnew3=[xnew2;xnew2];

figure(1);
plot(xnew1)
title('filt');

Fs = 122;
t = 0:1/Fs:1-1/Fs;
[P1,Q1] = rat(250/Fs);
abs(P1/Q1*Fs-250)

xnew2 = resample(xnew1,P1,Q1)

Fs1 = 250;
t = 0:1/Fs1:1-1/Fs1;
[P,Q] = rat(500/Fs1);
abs(P/Q*Fs1-500)

xnew3 = resample(xnew2,P,Q)

Fs1 = 500;
t = 0:1/Fs1:1-1/Fs1;
[P2,Q2] = rat(1000/Fs1);
abs(P/Q*Fs1-1000)

xnew4 = resample(xnew3,P2,Q2)


fs = 1000;
y2_1 = medfilt1(xnew4,fs/5);
y2_2 = medfilt1(y2_1,(3*fs)/5);
y2 = xnew4-y2_2;
title('median');

fid = fopen('PiiTech_mitdb_109_14_19.txt','w');
 fprintf(fid,'%d\n',y2);
 fclose(fid); 

figure 
plot(y2)
 
  close all
ecg_signal=y2(1:end);
figure(6)
plot(ecg_signal);


hold on;
close all 
% for i=1:1500;
     if (ecg_signal(1:1300)<0.08)
            ecg_signal=ecg_signal(1300:end)
     end  
figure(6)
plot(ecg_signal);
%%%%%%%%%%%%%%%% The below code is to get the DWT coefficients of the signal ------------- START %%%%
%------------ Resolution Level - 1 ------------------
[cA_l1, cD_l1] = dwt(ecg_signal,'haar');

%------------ Resolution Level - 2 ------------------
[cA_l2, cD_l2] = dwt(cA_l1,'haar');

%---Resolution level - 3 --------------------
[cA_l3, cD_l3] = dwt(cA_l2,'haar');

%---------- Resolution level - 4 ------------------
[cA_l4, cD_l4] = dwt(cA_l3,'haar');

%---------- Resolution Level - 5 -----------------
[cA_l5, cD_l5] = dwt(cA_l4,'haar');
%%%%%%%%%%%%%%%% The below code is to get the DWT coefficients of the signal ------------- END %%%%

%%%%%%%%%%%% Finding the threshold ____________________________________________________________ START____%%%%%
m=1;
n=175;
max_chunk=[];
for i=1:(length(cD_l3)/175)
     max_chunk=[max_chunk;max(cD_l3(m:n))];
     m=n+1;
n=n+175;
end
min_chunk=min(max_chunk)
T=min_chunk*0.45;

%%%%%%%%%%%% Finding the threshold ____________________________________________________________ END____%%%%%

cD_l3_logic=[];

%%%%%% Storing logic 0 and 1 in the memory_______________________________________________________Start_____%%%%%
for i=1:length(cD_l3)
    if(abs(cD_l3(i))>=T)
       cD_l3_logic(i)=1;
   else
    cD_l3_logic(i)=0;
   end
end
cD_l3_logic(i)=cD_l3_logic(i);
 
 for i=1:length(cD_l3)
     if(i<length(cD_l3)-2)
    if(cD_l3_logic(i)==1)
           if(cD_l3_logic(i)==1&cD_l3_logic(i+1)==1)
              cD_l3_logic(i)=0;
           
          else 
               cD_l3_logic(i)=cD_l3_logic(i); 
        end
      else
    cD_l3_logic(i)=cD_l3_logic(i);
    end
 end
 end
%%%%%% Storing logic 0 and 1 in the memory_______________________________________________________ END_____%%%%%
  j=1;
  k=1;
  
%%%%%% Finding R-peak index_______________________________________________________________________Start_____%%%%%
for i=1:length(cD_l3)
   if(cD_l3_logic(i)==1 && i<(length(cD_l3)-10) && i>10)
      
        max_index(j)=i;
       [value(k),min_index(k)]=min(cD_l3((max_index(j) - 10) : (max_index(j)+10)));
        if(min_index(k)>=10)
           min_index(k)= max_index(j)+(min_index(k)-10);
        else 
           min_index(k)= max_index(j)-(11-min_index(k));
        end
        
        if(max_index(j)<min_index(k))
          [R_peak_value(j),R_peak_index(j)]= min(ecg_signal((max_index(j)*8):(min_index(k)*8)));
               R_peak_index(j)=R_peak_index(j)+(max_index(j)*8)-1;
        else 
            [R_peak_value(j),R_peak_index(j)]= max(ecg_signal((min_index(j)*8):(max_index(k)*8)));
               R_peak_index(j)=R_peak_index(j)+(min_index(j)*8)-1;
       end
        
       j=j+1;
       k=k+1;
       
   end
end
R_peak_index=unique(R_peak_index,'first');

%%%%%% Finding R-peak index_________________________________________________________________________ End_____%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% eleminating the R-peaks whose distance are less than 290 Samples_____________________________START__________
i=1;
while i<=(length(R_peak_index)-1)
    if ((R_peak_index(i+1)- R_peak_index(i))<=290)
        if (abs(ecg_signal(R_peak_index(i)))>(0.02*abs(ecg_signal(R_peak_index(i+1)))))
               R_peak_index(i+1)=0;
                        i=i+2;
        elseif (abs(ecg_signal(R_peak_index(i)))<(0.02*abs(ecg_signal(R_peak_index(i+1)))))
            R_peak_index(i)=0;
             i=i+2;
        end

    else 
        R_peak_index(i)=R_peak_index(i)
        i=i+1;
    end
end 
R_peak_index(find(R_peak_index==0)) = [];
%%%%%%% eleminating the R-peaks whose distance are less than 250 Samples ______________________________END__________


i=1;
%%%%%%% eleminating the R-peaks whose distance are less than 250 Samples_____________________________START__________
%% The below code with the condition of 0.02 is changed for the data "mitdb_203_1788_1793"
while i<(length(R_peak_index)-1)
    if ((R_peak_index(i+1)- R_peak_index(i))<=290) 
        if (abs(ecg_signal(R_peak_index(i)))>(0.02*abs(ecg_signal(R_peak_index(i+1)))))
               R_peak_index(i+1)=0;
                        i=i+2;
        elseif (abs(ecg_signal(R_peak_index(i)))<(0.02*abs(ecg_signal(R_peak_index(i+1)))))
            R_peak_index(i)=0;
             i=i+2;
        end
    else 
        R_peak_index(i)=R_peak_index(i);
        i=i+1;
    end
end 
R_peak_index(find(R_peak_index==0)) = [];
%%%
if (R_peak_index(end)-R_peak_index(end-1)<=290)
    if (abs(ecg_signal(R_peak_index(end-1)))>(0.02*abs(ecg_signal(R_peak_index(end)))))
               R_peak_index(end)=0;
    end
end
%%%%%%% eleminating the R-peaks whose distance are less than 250 Samples_____________________________END__________


i=1;
%%%%%%% eleminating the R-peaks whose distance are less than 250 Samples_____________________________START__________
while i<(length(R_peak_index)-1)
    if ((R_peak_index(i+1)- R_peak_index(i))<290)
        R_peak_index(i+1)=0;
             i=i+2;
    else 
        R_peak_index(i)=R_peak_index(i)
        i=i+1;
    end
end 
R_peak_index(find(R_peak_index==0)) = [];
%%%%%%% eleminating the R-peaks whose distance are less than 250 Samples_____________________________END__________

%%%%%%%%%% finding the correct Max value for R-peak in the range of -56 to +56 ____________ Start____________
for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end

%%%%%%%%%% finding the correct Max value for R-peak in the range of -56 to +56 ____________ END____________
close all

%%% IF the R-peak to R-peak is >1200, then there might be few chances that some R-peaks may miss---So, take the range as mentioned in for loop shown and find 
%% find the 50% of max of R-peak and add that index in the R-peak array______________ START________
s=0
z=1;
for i=1:length(R_peak_index)-1
    if(R_peak_index(i+1)- R_peak_index(i))>1200
        for k= (R_peak_index(i)+50): (R_peak_index(i+1)-50)
            if (ecg_signal(k)>= (0.37*max(abs(ecg_signal(R_peak_index)))))%% 0.37 is added for the database "mitdb_223_1148_1153"
                       index9(1)=R_peak_index(i);
                      index9(z+1) =k;
                      z=z+1;
                      s=1;
                      break
            end
        end 
    end
    
end
%%% IF the R-peak to R-peak is >1200, then there might be few chances that some R-peaks may miss---So, take the range as mentioned in for loop shown and find 
%% find the 50% of max of R-peak and add that index in the R-peak array______________ END________
if(s)
index8=index9(2:end);
   [h j]= max(ecg_signal(index8));
   j=index8(j);
   
   if(j-index9(1)>400)
       R_peak_index=[R_peak_index, j];
       R_peak_index= sort(R_peak_index)
       
   end
end

s=0
z=1;
for i=1:length(R_peak_index)-1
    if(R_peak_index(i+1)- R_peak_index(i))>1200
        for kz= (R_peak_index(i)+50): (R_peak_index(i+1)-50)
            if (abs(ecg_signal(kz))>= (0.37*max(abs(ecg_signal(R_peak_index)))))%% 0.37 is added for the database "mitdb_223_1148_1153"
                       index9z(1)=R_peak_index(i);
                      index9z(z+1) =kz;
                      z=z+1;
                      s=1;
                      break
            end
        end 
    end
    
end
%%% IF the R-peak to R-peak is >1200, then there might be few chances that some R-peaks may miss---So, take the range as mentioned in for loop shown and find 
%% find the 50% of max of R-peak and add that index in the R-peak array______________ END________
if(s)
index8z=index9z(2:end);
   [hz jz]= max(ecg_signal(index8z));
   jz=index8z(jz);
   
   if(jz-index9z(1)>400)
       R_peak_index=[R_peak_index, jz];
       R_peak_index= sort(R_peak_index)
       
   end
end

%%%%%%%%%%%%%%%%%% Not getting why this logic has been written -- try giving an ECG signal having > 1200 HRV and analyse 
%how the arrays values are getting %changed_____________________________END____________
  
%%%%%%%% The below code is to find any missing R-peaks between the two R-Peaks, here we check the condition as shown and 
%%take -50 and +50 and find the 0.5 of max of R-peak in ecg. ______________ START ___________
zm=1;
sm=0;
for i=1:length(R_peak_index)-2
    if (R_peak_index(i+2)-R_peak_index(i+1))>= (1.8*(R_peak_index(i+1)-R_peak_index(i)))
         for k= (R_peak_index(i+1)+50):(R_peak_index(i+2)-50)
            if (ecg_signal(k)>= (0.5*max(abs((ecg_signal(R_peak_index))))))
                       index10(1)=R_peak_index(i+1);
                      index10(zm+1) =k;
                      zm=zm+1;
                      sm=1;
            end
        end 
    end
    
end
%%%%%%%% The below code is to find any missing R-peaks between the two R-Peaks, here we check the condition as shown and 
%%take -50 and +50 and find the 0.5 of max of R-peak in ecg. ______________ E ___________

if(sm)
index11=index10(2:end);
   [h j1]= max(ecg_signal(index11));
   j1=index11(j1);
   
   if(j1-index10(1)>100)
      R_peak_index=[R_peak_index, j1];
       R_peak_index= sort(R_peak_index)
       
   end
end

i=1;
while i<=(length(R_peak_index)-1)
    if ((R_peak_index(i+1)- R_peak_index(i))<250)
        R_peak_index(i+1)=0;
        if(i<length(R_peak_index)-2)
             i=i+2;
        else 
            i=i+1;
        end
    elseif(i>=length(R_peak_index)-2)
             i=i+1;
    else 
        R_peak_index(i)=R_peak_index(i)
        i=i+1;
    
    end
end 
R_peak_index(find(R_peak_index==0)) = [];



for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end


 
for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end



%%%%%%%%%%% Code changed--- Modified%%%%%%
% for i=1:length(R_peak_index)-1
%     if (abs(ecg_signal(R_peak_index(i+1))) <  (0.2 * abs(ecg_signal(R_peak_index(i)))))
%             R_peak_index(i+1)=0;
%     end
% end
% R_peak_index(find(R_peak_index==0)) = [];

% i=1;
while i<(length(R_peak_index)-1)
    if (abs(ecg_signal(R_peak_index(i+1))) <  (0.2 * abs(ecg_signal(R_peak_index(i)))))
        R_peak_index(i+1)=0;
             i=i+2;
    end
end
 R_peak_index(find(R_peak_index==0)) = [];


%%%%%%%%%%% Code changed--- Modified%%%%%%

lm=0;
ii=1;
index9=0;
km=1;
for ii=1: length(R_peak_index)-2
    if (abs(ecg_signal(R_peak_index(ii+1))) <  (0.67 * abs(ecg_signal(R_peak_index(ii)))))
                   [value9 index9] =   max(abs(ecg_signal((R_peak_index(ii+1)-10):(R_peak_index(ii+1)+150))))
                   if(index9>10)
                      index9= R_peak_index(ii+1)+index9-10;
                   else 
                       index9= R_peak_index(ii+1)-index9;
                   end
                 index12(km)= index9;
                 km=km+1;
                 R_peak_index=[R_peak_index(1:ii),index9,R_peak_index(ii+2:end)]
                   lm=1;
    end
end

if (lm)
 R_peak_index=sort(R_peak_index);
end

for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end



% lmn=0;
% iii=1;
% index99=0;
% kmn=1;
% for iii=1: length(R_peak_index)-2
% %     if (abs(ecg_signal(R_peak_index(iii+1))) <  (0.67 * abs(ecg_signal(R_peak_index(ii)))))
%                    [value99 index99] =   max(abs(ecg_signal((R_peak_index(iii+1)-10):(R_peak_index(iii+1)+196))))
%                    if(index99>10)
%                       index99= R_peak_index(iii+1)+index99-10;
%                    else 
%                        index99= R_peak_index(iii+1)-index99;
%                    end
%                  index19(km)= index99;
%                  kmn=kmn+1;
%                  R_peak_index=[R_peak_index(1:iii),index99,R_peak_index(iii+2:end)]
%                    lmn=1;
%     
% end
% 
% if (lm)
%  R_peak_index=sort(R_peak_index);
% end

for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end



for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_peak_index_dummy=[];


for i=1:length(R_peak_index)
    if (i<length(R_peak_index)) && (ecg_signal(R_peak_index(i))>0) && (R_peak_index(i+1)-R_peak_index(i)>700)
    
%     for xz=R_peak_index(i):R_peak_index(i)+180
        if(i<length(R_peak_index))
          R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):(R_peak_index(i)+200)))>=ecg_signal(R_peak_index(i)),1,'last');
          R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
          
%         else 
%             if(length(ecg_signal)-R_peak_index(end)>=180) && (ecg_signal(R_peak_index(end))>0)
%                 R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):(R_peak_index(i)+180)))>=ecg_signal(R_peak_index(i)),1,'last');
%                 R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
%             else   
%                 R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):length(ecg_signal)))>=ecg_signal(R_peak_index(i)),1,'last');
%                 R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
%            end
        end
    elseif (i==length(R_peak_index)) 
        if(length(ecg_signal)-R_peak_index(end)>=200) && (ecg_signal(R_peak_index(end))>0)
                R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):(R_peak_index(i)+200)))>=ecg_signal(R_peak_index(i)),1,'last');
                R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
             else   
                 R_peak_index_dummy(i) =R_peak_index(i);
%                R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
              end
        
        
    else 
       R_peak_index_dummy(i)=R_peak_index(i);
    end
end

R_peak_index=R_peak_index_dummy;






for i=1:length(R_peak_index)
    if(R_peak_index(1)>56) && (i==1)
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
            
    elseif (i==length(R_peak_index) && (5000-R_peak_index(end)>=56))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
    elseif (i>1) && (i<length(R_peak_index))
        
        [value,index]=max(abs(ecg_signal(R_peak_index(i)-56:R_peak_index(i)+56)));
          if (index>=56)
            index=index+R_peak_index(i)-56;
          else 
             index=R_peak_index(i)-(56-index)
          end
            R_peak_index(i) =index;
            R_Peak_value(i)=ecg_signal(index);
end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% R_peak_index_dummy=[];
% 
% 
% for i=1:length(R_peak_index)
%     if (i<length(R_peak_index)-1) && (ecg_signal(R_peak_index(i))>0) && (R_peak_index(i+1)-R_peak_index(i)>700)
%     
% %     for xz=R_peak_index(i):R_peak_index(i)+180
%         if(i<length(R_peak_index))
%           R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):(R_peak_index(i)+200)))>=ecg_signal(R_peak_index(i)),1,'last');
%           R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
%           
%         else 
%             if(length(ecg_signal)-R_peak_index(end)>=180) && (ecg_signal(R_peak_index(end))>0)
%                 R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):(R_peak_index(i)+200)))>=ecg_signal(R_peak_index(i)),1,'last');
%                 R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
%             else   
%                 R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):length(ecg_signal)))>=ecg_signal(R_peak_index(i)),1,'last');
%                 R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
% end
% end
%     elseif (i==length(R_peak_index)) 
%         if(length(ecg_signal)-R_peak_index(end)>=180) && (ecg_signal(R_peak_index(end))>0)
%                 R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):(R_peak_index(i)+200)))>=ecg_signal(R_peak_index(i)),1,'last');
%                 R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
%             else   
%                 R_peak_index_dummy(i) =find(abs(ecg_signal(R_peak_index(i):length(ecg_signal)))>=ecg_signal(R_peak_index(i)),1,'last');
%                 R_peak_index_dummy(i)=R_peak_index_dummy(i)+ R_peak_index(i);
% end
%         
%         
%     else 
%         R_peak_index_dummy(i)=R_peak_index(i);
%         end
% end
% 
R_peak_index_dummy=R_peak_index;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%% R-peak is always +ve (positive), the below code makes it%%
 for ii=1: length(R_peak_index)
    if(ecg_signal(R_peak_index(ii))<0 && ii==1)
        if(R_peak_index (1)>89)%%%%%% the value 79 is working for all the signals except "mitdb_217_181_186" and "mitdb_217_183_188", so the value is chnaged to "88"
%        [value99,index99] =max(ecg_signal(R_peak_index(ii)-79):ecg_signal(R_peak_index(ii)));
            [value99,index99]=max((ecg_signal(R_peak_index(ii)-89:R_peak_index(ii))));
            R_peak_index(ii)=(R_peak_index(ii)-89)+index99
        else 
            [value99,index99]=max((ecg_signal(1:R_peak_index(ii))));
            R_peak_index(ii)=index99
        end
    elseif (ecg_signal(R_peak_index(ii))<0)
        [value99,index99]=max((ecg_signal(R_peak_index(ii)-89:R_peak_index(ii))));
        R_peak_index(ii)=(R_peak_index(ii)-89)+index99
    end
 end
 
 
for i=1:length(R_peak_index)
    R_Peak_value99(i)=ecg_signal(R_peak_index(i));
end


for i=1 : length(R_peak_index)-1
    HRV(i)= (R_peak_index(i+1)- R_peak_index(i));
    
end



%%%%%%%%%%%%%%%% Boundary Code Started   %%%%%%%%%%%%%%
%% The command "((((R_peak_index(i+1)- R_peak_index(i))/2) + R_peak_index(i))+190<R_peak_index(i+1))" is used to avoid two R-peaks in a single boundary.... This 
%% usually happens when the R-peak to R-peak distance is very less... i.e in the database "mitdb_205_296_301" the distance between the 7 th and 6 th beat is 
% very less---- so, the mid point of R7 and R6 + 190 will exceed the R7
% value...to avoid this we have written the above command
%  for i=1:length(R_peak_index)-1
%      %I=5
%     if ((R_peak_index(i+1)- R_peak_index(i))<650) && ((((R_peak_index(i+1)- R_peak_index(i))/2) + R_peak_index(i))+190<R_peak_index(i+1))
%          for k=1:(round(((R_peak_index(i)+R_peak_index(i+1))/2)):(round(((R_peak_index(i)+R_peak_index(i+1))/2))+190));
%                 if ((ecg_signal(k+round((R_peak_index(i)+R_peak_index(i+1))/2)) * ecg_signal(k+1+round((R_peak_index(i)+R_peak_index(i+1))/2)))<0) || (ecg_signal(k+round((R_peak_index(i)+R_peak_index(i+1))/2))==0)
%                        boundary(i)= k + (round((R_peak_index(i)+R_peak_index(i+1))/2));
%                        break;
%                 else 
%                  end
%          end
%     elseif (R_peak_index(i+1)- R_peak_index(i)<190)
%             boundary(i)=((R_peak_index(i+1)+R_peak_index(i))/2);
%     else
%         boundary(i)=((R_peak_index(i+1)+R_peak_index(i))/2)+40;
%         
%      end
%  end 
 
 


for i=1:length(R_peak_index)-1
     %I=5
    if ((R_peak_index(i+1)- R_peak_index(i))<658) && ((((R_peak_index(i+1)- R_peak_index(i))/2) + R_peak_index(i))+190<R_peak_index(i+1))
         for k=round(((R_peak_index(i)+R_peak_index(i+1))/2)):(round(((R_peak_index(i)+R_peak_index(i+1))/2))+190);
                if ((ecg_signal(k) * ecg_signal(k+1))<0) || (ecg_signal(k)==0)
                       boundary(i)= k ;
                       break;
                else 
                    boundary(i)=((R_peak_index(i+1)+R_peak_index(i))/2);
                 end
         end
    elseif (R_peak_index(i+1)- R_peak_index(i)<190)
            boundary(i)=((R_peak_index(i+1)+R_peak_index(i))/2);
    else
        boundary(i)=((R_peak_index(i+1)+R_peak_index(i))/2)+40;
        
     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:length(HRV)-1
 if (abs(HRV(i+1)-HRV(i))>=100) && ((((R_peak_index(i+1)- R_peak_index(i))/2) + R_peak_index(i))+190<R_peak_index(i+1));
% if (abs(R_peak_index(i+1)-R_peak_index(i))>=100);
    
         for k= round(((R_peak_index(i)+R_peak_index(i+1))/2)): (round(((R_peak_index(i)+R_peak_index(i+1))/2))+190);
                if ((ecg_signal(k) * ecg_signal(k+1)<0) || (ecg_signal(k)==0))
                       boundary(i)= k ;
                       break;
                 end
           
        end
    else
        boundary(i)=boundary(i);
        
     end
 end 
 
 boundary=round(boundary);
  
 if (R_peak_index(1)>= (R_peak_index(2)- boundary(1)) &&  R_peak_index(2)-R_peak_index(1)>500)
    Boundary_start=(R_peak_index(2)- boundary(1))
    Boundary_start=R_peak_index(1)-Boundary_start;
    gridxy(Boundary_start,'Color','red','linewidth',2);
     boundary = [Boundary_start boundary];
end

if((length(ecg_signal)-R_peak_index(end))>= (boundary(end)-R_peak_index(end-1)))
    boundary_end=((boundary(end)-R_peak_index(end-1)));
    boundary_end=boundary_end+R_peak_index(end);
    gridxy(boundary_end,'Color','m','linewidth',2);
    boundary = [boundary boundary_end]
 end

 boundary=round(boundary)
 
 
  
 
 
 %%%%%% setting the proper boundaries if the boundary index value is >0.03
 for kz=1:length(boundary)
    if(ecg_signal(boundary(kz))>0.03 && kz<length(boundary))    
        for kk= boundary(kz):boundary(kz)+74
            if (kk<length(ecg_signal)) && ((ecg_signal(kk)* ecg_signal(kk+1))<0 || ecg_signal(kk)<0.015)
                        boundary(kz)=kk;
            end
        end
    elseif kz==length(boundary) && (length(ecg_signal)-boundary(end)<=100) && abs((ecg_signal(boundary(kz))>0.03))
         for kk= boundary(kz):length(ecg_signal)
            if (kk<length(ecg_signal)) && ((ecg_signal(kk)* ecg_signal(kk+1))<0 || ecg_signal(kk)<0.02 || ecg_signal(kk+1)<0.02)
                        boundary(kz)=kk;
            end
         end
     elseif kz==length(boundary) && (length(ecg_signal)-boundary(end)>100) && abs(ecg_signal(boundary(kz)))>0.03 && ecg_signal(boundary(end)+1)> ecg_signal(boundary(end))
         for kk= round((R_peak_index(end)+boundary(end))/2):boundary(end)+80;
            if (kk<length(ecg_signal)) && ((ecg_signal(kk)* ecg_signal(kk+1))<0 || ecg_signal(kk)<0.015)
                        boundary(kz)=kk;
            end
         end
    end
 end

 
 if (ecg_signal(boundary(end)>1))
     if  boundary(end)<length(ecg_signal) && ecg_signal(boundary(end)+1)>ecg_signal(boundary(end)) && (boundary(end)-R_peak_index(end)>400)
         for  kk= boundary(end)-125:boundary(end)
             if (kk<length(ecg_signal)) && ((ecg_signal(kk)* ecg_signal(kk+1))<0 || ecg_signal(kk)<0.015)
                        boundary(end)=kk;
             end
         end
     end
     
 end
     
     
 for i=1:length(boundary)
     boundary_value(i)=ecg_signal(boundary(i));
 end
 
 
  figure
  plot(ecg_signal)

clc
close all
boundary=boundary'
R_peak_index=R_peak_index'
HRV=HRV'
figure
plot(ecg_signal)
hold on
plot(R_peak_index,R_Peak_value99,'rx','linewidth',10);
plot(boundary,boundary_value,'cx','linewidth',10);


% s=strcat(name,'-',num2str(i));
%     fullFileName = fullfile(path,s);
%     saveas(gcf, fullFileName);
%     close(gcf);
    
