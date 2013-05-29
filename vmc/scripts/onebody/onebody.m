function [count centers] = onebody(filename, numprocs, nParticles, nSamples, binary)
    % loading data from all processors into "data"
    data = [];
    endings = "xyzr";
    input = [];
    data = [];
    if (binary)
        for i=0:(numprocs-1)
            xyzr = [];
            for j=1:4
                fileName = [filename '_' num2str(i) "of" num2str(numprocs)  "_nParticles" num2str(nParticles) "_nSamples" num2str(floor(nSamples/numprocs), 10) "_" endings(j) ".bin"];
                if (exist([fileName], 'file') ~= 2)
                    disp(['Wrong filename, "', fileName, '" does not exist!']);
                    return
                end
                [fileID message] = fopen(fileName);
                columnIn = fread(fileID, [floor(nSamples/numprocs)*nParticles, 1], 'double=>double');
                xyzr = [xyzr columnIn];
                
                %size(columnIn)
                %size(xyzr)
                %floor(nSamples/numprocs)*nParticles
            end
            data = [data; xyzr];
            %size(data)
        end
    else
        for i=0:(numprocs-1)
            fileName = [filename '_' num2str(i) "of" num2str(numprocs)  "_nParticles" num2str(nParticles) "_nSamples" num2str(nSamples/numprocs) ".dat"];
            input = load(fileName);
        end
		data = [data; input];
	end
    clear input;
    clear xyzr;
    clear columnIn;
    
	[a b] = size(data);
    N = a;
    
    size(data)
    
    figure;
    [count centers] = hist(data(:,4), 500);
    fac = (centers(2) - centers(1)).*centers.*centers;
    plot(centers, count./fac, '.')
    
    %trapz(centers, count)
    
    figure;
    plot(centers, count, '.')
    
    %figure;
    %[count centers] = hist(data(:,1), 500);
    %plot(centers, count, '.')
    %%%%plot(centers, count.*(centers.*centers), '.')

