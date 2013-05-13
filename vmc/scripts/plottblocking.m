function plottblocking(filename)
	input = load(filename);
	data = input.results;
	blockSize = data(:,1);
	sigma = data(:,3);
	
	figure;
	plot(blockSize, sigma);
	xlabel('BlockSize');
	ylabel('\sigma', 'interpreter', 'tex');
	
