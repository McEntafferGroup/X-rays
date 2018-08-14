function [ imagex ] = add_xrays( image, num, peak, sigma, show )
% image: the gaussian image to add x-rays to
% num:   the number of events to add
% peak:  counts in adu of x-rays
% sigma: the sigma of the x-ray events
% show:  boolean of whether to display the image or not    
%     
% returns: 2D array of image pixels

energy = random('norm',peak,sigma,num,1);

[ Y, X ] = size(image);

x_pos  = randi(X,num,1);
y_pos  = randi(Y,num,1);

Single_Events = table(x_pos,y_pos,energy)

imagex = image;

for i = 1:num
    imagex(y_pos(i),x_pos(i)) = imagex(y_pos(i),x_pos(i)) + energy(i);
end

if (show)
    figure;
    colormap('hot');
    imagesc(imagex);
    colorbar;
end

end