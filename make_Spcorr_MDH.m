function make_Spcorr_MDH(L,Jstr,Jdis,Jz,Pdist,Jseed)
%make_Spcorr_MDH(L,Jstr,Jdis,Jz,Pdist,Jseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_Spcorr_MDH
% function to make spin corr from components
% S.S = SzSz + 0.5*(SpSm + SmSp)
% 
% Andrew Goldsborough - 27/01/2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('printing Spcorr\n');

%open files to read in data
fnameSz = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_Szcorr_MDH.txt');
Spcorr = importdata(fnameSz);

fnameSpSm = strcat('./SpSmcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_SpSmcorr_MDH.txt');
corr = importdata(fnameSpSm);
Spcorr(:,3) = Spcorr(:,3) + 0.5*corr(:,3);

fnameSmSp = strcat('./SmSpcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_SmSpcorr_MDH.txt');
corr = importdata(fnameSmSp);
Spcorr(:,3) = Spcorr(:,3) + 0.5*corr(:,3);

%open files to write to
fnameSp = strcat('./Spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_Spcorr_MDH.txt');
fprintf(strcat(fnameSp,'\n'));
fidSpcorr = fopen(fnameSp, 'w');

%print to file
for i=1:size(Spcorr,1)
    fprintf(fidSpcorr,'%d %d %.15e\n',Spcorr(i,1),Spcorr(i,2),Spcorr(i,3));
end

%close file
fclose(fidSpcorr);