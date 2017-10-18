function MDH_PBC(L,Jstr,Jdis,Jz,Pdist,Jseed)
%MDH_PBC(L,Jstr,Jdis,Jz,Pdist,Jseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MDH PBC
% build MDH wavefunction and calculate corr and ee
% based on 10.1103/PhysRevB.22.1305 by Dasgupta and Ma
% currently only for XX and XXX
%
% Andrew Goldsborough - 26/01/2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%inputs
% L = 10;         %chain length
% Jstr = 1;       %overall J strength
% Jdis = 2;       %disorder strength
% Jz = 1;         %anisotropy
% Pdist = 6;      %coupling distribution
% Jseed = 2;      %seed for rng, 0 => shuffle

%coupling distribution
%0 => manual
%1 => 2 theta(K-1/2)
%2 => 1
%3 => uniform around Jstr normalised by Jstr
%4 => uniform around Jstr un-normalised Jdis
%5 => box distribution of Hikihara AF (10.1103/PhysRevB.60.12116)
%6 => Laflorencie's infinite disorder distribution (10.1103/PhysRevB.72.140408)

if (Jz ~= 0) && (Jz ~= 1)
    error('currently only for XX and XXX models');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate couplings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set seed, set by clock if zero
if Jseed > 0
    rng(Jseed);
else
    rng('shuffle');
end

%set the probability distribution
if Pdist==0
    %custom set
    J = [0.1,Jdis,0.1,1,Jdis,1,0.1,Jdis,0.1,0.01];
elseif Pdist==1
    %P(K) = 2 theta(K-1/2)
    J = zeros(1,L) + Jstr*random('unif',0.5,1,[1,L]);
elseif Pdist==2
    %P(K) = 1
    J = zeros(1,L) + Jstr*(rand(L,1));
elseif Pdist==3
    %uniform around Jstr normalised by Jstr
    J = zeros(1,L) + Jstr + Jstr*Jdis*(rand(1,L) - 0.5);
elseif Pdist==4
    %uniform around Jstr un-normalised Jdis
    J = zeros(1,L) + Jstr + Jdis*(rand(1,L) - 0.5);
elseif Pdist==5
    %box distribution of Hikihara AF
    J = zeros(1,L) + Jdis*random('unif',0,1,[1,L]);
elseif Pdist==6
    %Laflorencie's infinite disorder distribution
    J = rand(1,L).^Jdis;
end

rng('shuffle');

%print inputs
fprintf('L = %d, Jstr = %f, Jdis = %f, Jz = %f, Pdist = %d, Jseed = %d\n',L,Jstr,Jdis,Jz,Pdist,Jseed);

%print J to file
fprintf('printing interaction strengths\n');
fname = strcat('./J/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_J.txt');
fidJ = fopen(fname, 'w');
for i = 1:L
    fprintf(fidJ,'%.15e\n',J(i));
end
fclose(fidJ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do MDH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energy = 0;
sites = 1:L;
connections = zeros(1,L); %sites connected by singlet

for iteration = 1:(L/2)-1
    
    %find the maximum J
    [Jmax,Jmaxpos] = max(abs(J));
    
    %save site numbers connected by singlet
    connections(sites(PBC_pos(Jmaxpos,size(J,2)))) = sites(PBC_pos(Jmaxpos+1,size(J,2)));
    connections(sites(PBC_pos(Jmaxpos+1,size(J,2)))) = sites(PBC_pos(Jmaxpos,size(J,2)));
    
    %calculate the singlet energy
    if Jz == 1
        energy = energy - 0.75*Jmax - (3/(16*Jmax))*(J(PBC_pos(Jmaxpos-1,size(J,2)))^2 + J(PBC_pos(Jmaxpos+1,size(J,2)))^2);
    elseif Jz == 0
        energy = energy - 0.5*Jmax;
    end
    
    %renormalise coupling
    J(PBC_pos(Jmaxpos-1,size(J,2))) = J(PBC_pos(Jmaxpos-1,size(J,2)))*J(PBC_pos(Jmaxpos+1,size(J,2)))/((1+Jz)*J(Jmaxpos));
    
    %remove the two sites
    J(Jmaxpos) = [];
    J(PBC_pos(Jmaxpos,size(J,2))) = [];
    sites(Jmaxpos) = [];
    sites(PBC_pos(Jmaxpos,size(sites,2))) = [];
end

%final iteration
Jmaxpos = 1;
connections(sites(PBC_pos(Jmaxpos,size(J,2)))) = sites(PBC_pos(Jmaxpos+1,size(J,2)));
connections(sites(PBC_pos(Jmaxpos+1,size(J,2)))) = sites(PBC_pos(Jmaxpos,size(J,2)));

if Jz == 1
    energy = energy - 0.75*(J(1)+J(2));
elseif Jz == 0
    energy = energy - 0.5*(J(1)+J(2));
end

%print energy to file
fprintf('printing energy\n');
fname = strcat('./energy/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_energy_MDH.txt');
fidenergy = fopen(fname, 'w');
fprintf(fidenergy,'%.15e\n',energy);
fclose(fidenergy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open files to print to
fprintf('printing correlation functions\n');
fname = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_Szcorr_MDH.txt');
fidSzcorr = fopen(fname, 'w');

fname = strcat('./SpSmcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_SpSmcorr_MDH.txt');
fidSpSmcorr = fopen(fname, 'w');

fname = strcat('./SmSpcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_SmSpcorr_MDH.txt');
fidSmSpcorr = fopen(fname, 'w');

for i = 1:L-1
    for j = i+1:L
        if connections(i) == j
            %Szcorr = -0.25 if directly connected
            fprintf(fidSzcorr,'%d %d %.15e\n',i,j,-0.25);
            
            %SpSmcorr = -0.5 if directly connected
            fprintf(fidSpSmcorr,'%d %d %.15e\n',i,j,-0.5);
            
            %SmSpcorr = -0.5 if directly connected
            fprintf(fidSmSpcorr,'%d %d %.15e\n',i,j,-0.5);
        else
            %zero otherwise
            fprintf(fidSzcorr,'%d %d %.15e\n',i,j,0);
            fprintf(fidSpSmcorr,'%d %d %.15e\n',i,j,0);
            fprintf(fidSmSpcorr,'%d %d %.15e\n',i,j,0);
        end
    end
end

fclose(fidSzcorr);
fclose(fidSpSmcorr);
fclose(fidSmSpcorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%entanglement entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open file to print to
fprintf('printing ee\n');
fname = strcat('./ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(Pdist),'_',num2str(Jseed),'_ee_MDH.txt');
fidee = fopen(fname, 'w');

%ee is the number of singlets that need to be broken to separate the chain
for i = 1:L-1
    
    %work from left to right
    bonds = 0;
    
    %shift connections to match i
    con_ee = circshift(PBC_pos(connections-i,L),[0 -i]);
    
    for j = i+1:L
        if con_ee(j-i) > (j-i)
            %singlet is being opened
            bonds = bonds + 1;
        else %con_ee(j-i) < (j-i)
            %singlet being closed
            bonds = bonds - 1;
        end
        
        fprintf(fidee,'%d %d %.15e\n',i,j,bonds);
    end
end

fclose(fidee);
toc
end
