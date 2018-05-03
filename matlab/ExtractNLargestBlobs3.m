% Function to return the specified number of largest or smallest blobs in a binary image.
% If numberToExtract > 0 it returns the numberToExtract largest blobs.
% If numberToExtract < 0 it returns the numberToExtract smallest blobs.
% Example: return a binary image with only the largest blob:
%   binaryImage = ExtractNLargestBlobs(binaryImage, 1)
% Example: return a binary image with the 3 smallest blobs:
%   binaryImage = ExtractNLargestBlobs(binaryImage, -3)
function mask = ExtractNLargestBlobs(mask, n)
try
	% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
	[labeledImage, numberOfBlobs] = bwlabeln(mask, 6);
	blobMeasurements = regionprops3(labeledImage, 'volume');
	% Get all the areas
	allVolumes = [blobMeasurements.Volume];
	if n > 0
		% For positive numbers, sort in order of largest to smallest.
		% Sort them.
		[sortedVolumes, sortIndexes] = sort(allVolumes, 'descend');
	elseif n < 0
		% For negative numbers, sort in order of smallest to largest.
		% Sort them.
		[sortedVolumes, sortIndexes] = sort(allVolumes, 'ascend');
		% Need to negate numberToExtract so we can use it in sortIndexes later.
		n = -n;
	else
		% numberToExtract = 0.  Shouldn't happen.  Return no blobs.
		mask = false(size(mask));
		return;
	end
	% Extract the "numberToExtract" largest blob(a)s using ismember().
	biggestBlob = ismember(labeledImage, sortIndexes(1:n));
	% Convert from integer labeled image into binary (logical) image.
	mask = biggestBlob > 0;
catch ME
	errorMessage = sprintf('Error in function ExtractNLargestBlobs().\n\nError Message:\n%s', ME.message);
	fprintf(1, '%s\n', errorMessage);
	%uiwait(warndlg(errorMessage));
end
