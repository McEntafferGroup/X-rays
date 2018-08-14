function [ image ] = event_consolidator( imagex, thresh )

[ Y, X ] = size(imagex);

bg   = mean2(imagex);
sigma = std2(imagex);

s = thresh*sigma;

image = (imagex - bg);

for i = 1:Y
    
    for j = 1:X
        
        pixel = image(i,j);
        
        if ( pixel > s )
            
            image = surrounding(image,i,j,s,sigma);
            
        end
        
    end
    
end

image = (image + bg);

end