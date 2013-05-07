function result = processblock(data, blockSize)
	samples = length(data);
	blocks = samples/blockSize;
	averageEnergyBlock = zeros(blocks,1);
	
	% computing the mean value of every block 
	for i=1:blocks
		% finding the average over one block
		%temp = data(((i-1)*blockSize+1):i*blockSize);
		averageEnergyBlock(i) = mean(data(((i-1)*blockSize+1):i*blockSize));
	end

	% calcualting the mean and variance of all the blocks
	E = mean(averageEnergyBlock);
	E2 = mean(averageEnergyBlock.*averageEnergyBlock);
	sigma = sqrt((E2 - E.*E)/blocks);
	
	result = zeros(1,3);
	result(1) = blockSize;
	result(2) = E;
	result(3) = sigma;
	
