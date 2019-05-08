close all;
fclose('all');

% read data
pth=input('Directory [default=current directory]  ');
	if isempty(pth)
        %pth=pwd;
        pth='H:\projects\research practicum\single molecule fluorescence\Matlab\HJA-data from Ivo'
   	end
cd(pth);

%cd('C:\user\tir data\May 25th 2004\Rep\Rep Hepes');

% read data
fname=input('index # of filename [default=1]  ');
	if isempty(fname)
   	%fname=1;
    fname=6;
	end
fname=num2str(fname);
['hel' fname '.traces']

% define time unit
%timeunit=0.1;
timeunit=1
%timeunit=input('Time Unit [default=0.1 sec]  ');
%	if isempty(timeunit)
%      timeunit=0.1;
%	end
 
maxcps=1500;
%maxcps=input('Maximun FRET value on the graph [default=2500]  ');
%	if isempty(maxcps)
%   	maxcps=2500;  %10 msec
%   end
   
leakage = 0;
%leakage=input('Donor leakage correction [default=0]  ');
%   if isempty(leakage)
%      leakage=0;
%   end
%   leakage   

fid=fopen(['hel' fname '.traces'],'r');

len=fread(fid,1,'int32');
disp('The len of the time traces is: ')
disp(len);
Ntraces=fread(fid,1,'int16');
%Ntraces=1298;
disp('The number of traces is:')
disp(Ntraces);


raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);
% convert into donor and acceptor traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
donor=zeros(Ntraces/2,len);
acceptor=zeros(Ntraces/2,len);
Data(index)=raw(index);
for i=1:(Ntraces/2),
   donor(i,:)=Data(i*2-1,:);
   acceptor(i,:)=Data(i*2,:);
end
   
time=(0:(len-1))*timeunit;
% calculate, plot and save average traces
dAvg=sum(donor,1)/Ntraces*2;
aAvg=sum(acceptor,1)/Ntraces*2;
avgOutput=[time' dAvg' aAvg'];
avgFileName=[fname '_avg.dat'];
save(avgFileName,'avgOutput','-ascii');
% calculate E level from the first 10 points and plot histograms of E level and total intensity. Also save the same info.
j=0;
elevel=zeros(Ntraces/2,1);
total=elevel;

output=[fname '_elevel_10p.dat'];
for i=1:(Ntraces/2);
   tempD=sum(donor(i,(3:12)),2);
   tempA=sum(acceptor(i,(3:12)),2);
   total(i)=(tempA+tempD)/10;
   elevel(i)=tempA/(tempA+tempD);
end


i=0;
while ((Ntraces/2 - 1) > 0) & (i < Ntraces/2 - 2) ,
   i = i + 1 ;
   % Trace window
   figure(1);
   subplot(2,1,1);

  plot(time,donor(i,:),'g',time,acceptor(i,:)-leakage*donor(i,:),'r');
   title(['  Molecule ' num2str(i)], 'fontsize',11,'fontweight','bold');
   ylabel ('Fluorescence Intensity (a.u.)', 'fontsize',10,'fontweight','bold');
   temp=axis;
   temp(3)=-maxcps/10;
   temp(4)=maxcps;
   grid on;
   axis(temp);
   zoom on;
   
   subplot(2,1,2);
   fretE=(acceptor(i,:)-leakage*donor(i,:))./(acceptor(i,:)+donor(i,:));
   for m=1:len,
      if acceptor(i,m)+donor(i,m)==0
         fretE(m)=-0.5;
      end
   end % This is to avoid undefined fretE interfering with future analysis
   
   plot(time,fretE,'b');
   xlabel ('Time (sec)','fontsize',10,'fontweight','bold');
   ylabel ('FRET efficiency (normalized)', 'fontsize',10,'fontweight','bold');
   temp=axis;
   temp(3)=-0.1;
   temp(4)=1.1;
   axis(temp);
   grid on;
   zoom on;
   
     
   ans=input('press s to save and press Enter to pass.','s');
   
   if ans=='g'
      nml = input('which molecule?');
      i=nml - 1;
       end

   if ans=='b'
      i=i - 2;
    end

   if ans=='s'
      j=j+1;
      mN(j)=i;
      fname1=[fname ' tr' num2str(i) '.dat'];
      output=[time' donor(i,:)' acceptor(i,:)'];
      save(fname1,'output','-ascii') ;
   end
   
   if ans=='x',
      fclose all;
	  clear all;
      break;
    end      
 
 end
