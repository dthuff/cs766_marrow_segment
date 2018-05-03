function [ output_args ] = plot_slice( scan , slice )

    figure();
    if max(scan(:)) == 1 %for masks
        imshow(squeeze(scan(:,:,slice)));
    else
        imshow(uint8(squeeze(scan(:,:,slice))));
    end        
end

