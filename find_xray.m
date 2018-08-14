function find_xray( image, thresh )
% image:  the image in which to find the events
% thresh: the threshold of sigms above which are events
%
% returns: void (will print a pretty table

peak  = mean2(image);
sigma =  std2(image);

sub = image - peak - (thresh*sigma);

[ y_pos, x_pos ] = find(sub > 0);

energy = zeros(length(y_pos),1);
for i = 1:numel(energy)
    energy(i) = image(y_pos(i), x_pos(i));
end

Found_Events = table(x_pos, y_pos, energy)

end