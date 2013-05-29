function blocking(filename, numprocs)
    deltaBlockSize = 10;
    maxBlockSizeThreshold = 1e4;
    stepLength = 10;
    outFilename = "blocking_results.mat";
    
    % loading data from all processors into "data"
    data = [];
	for i=0:(numprocs-1)
		fileName = [filename "_" num2str(i) "of" num2str(numprocs) ".dat"];
		%fileName = [filename "_" num2str(i) "of" num2str(numprocs) ".mat"];
		input = load(fileName);
		data = [data; input(:,1)];

		disp(['done loading process = ' num2str(i)]);
	end
	
	N = length(data);
	maxBlockSize = round(N/10);
	if (maxBlockSize > maxBlockSizeThreshold)
		maxBlockSize = maxBlockSizeThreshold;
	end

	disp('Done loading files');
	disp(['maxBlockSize/deltaBlockSize = ' num2str(maxBlockSize/deltaBlockSize)]);

	results = zeros(maxBlockSize/deltaBlockSize, 3);
	for i=1:maxBlockSize/deltaBlockSize
		blockSize = i*deltaBlockSize;
		results(i,:) = processblock(data, blockSize);
		disp(['i = ' num2str(i)]);
	end

	save(outFilename, 'results');

