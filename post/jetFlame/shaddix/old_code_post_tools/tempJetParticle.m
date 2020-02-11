%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guangyuan Sun 04/30/12
% Process mixing layer particle results
% Standalone code.
% Plots particle trajectories.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%%%%%%%%%%%%%%% get number of particles %%%%%%%%%%%%%%%%

command = ['../data/particles_init.dat'];
file = fopen(command);
ln = fgetl(file);
ln = fgetl(file);
k = 1;
while(~feof(file))
    ln = fgetl(file);
    save(k,:) = [sscanf(ln, '%f')]';
    k = k+1;
end
npart = length(save);

%%%%%%%%%%%%%%% get number of dump time steps %%%%%%%%%%%%%%%%

command = ['../input/dumpTimes.inp'];
file = fopen(command);
ln = fgetl(file);
k = 1;
while(~feof(file))
    ln = fgetl(file);
    dataTime(k,:) = [sscanf(ln, '%f')]';
    k = k+1;
end
ndumptime = length(dataTime);

for ifile = 1:ndumptime
    ifile

    command = ['../data/dmp_part_',num2str(ifile),'.dat']
    file = fopen(command);

    ln = fgetl(file);
    ln = fgetl(file);

    clear A;
    i = 1;
    while(~feof(file))
        ln = fgetl(file);
        A(i,:) = [sscanf(ln, '%f')]';
        i = i+1;
    end
    position(ifile,:) = A(:,2);
    fclose(file);
end

%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%

for i = 1:npart
    i
    hold on
    plot(position(:,i),dataTime(:),'k');
end
xlabel('ODT line position (mm)')
ylabel('Time (ms)')

