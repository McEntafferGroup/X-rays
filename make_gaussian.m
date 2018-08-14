function [ image ] = make_gaussian( size, peak, sigma, show )
% size:  (x,y)
% peak:  counts in adu
% sigma: the sigma of the gaussian noise
% show:  boolean of whether to display the image or not    
%     
% returns: 2D array of image pixels

x = size(1);
y = size(2);

image = random('norm',peak,sigma,y,x);

if (show)
    figure;
    colormap('hot');
    imagesc(image);
    colorbar;
    figure;
    histogram(image);
end

end