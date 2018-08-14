function [ imagex ] = add_split_events( image, num, type, peak, sigma, show )
% image: the image to add events to
% num:   number of events to be added
% type:  2,4,cross, or all splits (2,4,5,3 respectively)
% peak:  counts in adu of peak x-ray events
% sigma: the sigma of the x-ray events
% show:  boolean of whether to display the image or not
% 
% returns: 2D array of image pixels

imagex = image;

[ Y, X ] = size(image);
    

if (type == 2)
    
    xy = randi(2,1,num)-1;
    x = xy;
    y = ~x;
    
    energy = random('norm',peak,sigma,num,1);
    
    x_pos = zeros(num,1);
    y_pos = zeros(num,1);
    
    for i = 1:num
        
        x_pos(i) = randi(X-x(i));
        y_pos(i) = randi(Y-y(i));
        
        split1 = rand;
        split2 = 1-split1;
        
        imagex(y_pos(i),x_pos(i))           = imagex(y_pos(i),x_pos(i)) + (split1*energy(i));
        imagex(y_pos(i)+y(i),x_pos(i)+x(i)) = imagex(y_pos(i)+y(i),x_pos(i)+x(i)) + (split2*energy(i));
        
    end
    
    Two_Split_Event = table(x_pos,y_pos,energy)
    
end

if (type == 4)
    
    energy = random('norm',peak,sigma,num,1);
    
    x_pos = zeros(num,1);
    y_pos = zeros(num,1);
    for i = 1:num
        
        bit1 = random('norm',0.25,0.01);
        bit2 = random('norm',0.25,0.01);
        bit3 = random('norm',0.25,0.01);
        bit4 = random('norm',0.25,0.01);
        
        x_pos(i) = randi(X-1);
        y_pos(i) = randi(Y-1);
        
        imagex(y_pos(i),   x_pos(i))   = imagex(y_pos(i),   x_pos(i))   + bit1*energy(i);
        imagex(y_pos(i)+1, x_pos(i))   = imagex(y_pos(i)+1, x_pos(i))   + bit2*energy(i);
        imagex(y_pos(i),   x_pos(i)+1) = imagex(y_pos(i),   x_pos(i)+1) + bit3*energy(i);
        imagex(y_pos(i)+1, x_pos(i)+1) = imagex(y_pos(i)+1, x_pos(i)+1) + bit4*energy(i);
        
        % make it proper from all pixel values
        energy(i) = bit1*energy(i)+bit2*energy(i)+bit3*energy(i)+bit4*energy(i);
        
    end
    
    Four_Split_Event = table(x_pos,y_pos,energy)
    
end

if (type == 5)
    
    energy = random('norm',peak,sigma,num,1);
    
    bit0 = 3*random('norm',0.1429,0.01);
    bit1 = random('norm',0.1429,0.01);
    bit2 = random('norm',0.1429,0.01);
    bit3 = random('norm',0.1429,0.01);
    bit4 = random('norm',0.1429,0.01);
    
    x_pos = zeros(num,1);
    y_pos = zeros(num,1);
    
    for i = 1:num
        
        x_pos(i) = randi([2,X-1]);
        y_pos(i) = randi([2,Y-1]);
        
        imagex(y_pos(i),   x_pos(i))   = imagex(y_pos(i),   x_pos(i))   + bit0*energy(i);
        imagex(y_pos(i)+1, x_pos(i))   = imagex(y_pos(i)+1, x_pos(i))   + bit1*energy(i);
        imagex(y_pos(i),   x_pos(i)+1) = imagex(y_pos(i),   x_pos(i)+1) + bit2*energy(i);
        imagex(y_pos(i)-1, x_pos(i))   = imagex(y_pos(i)-1, x_pos(i)+1) + bit3*energy(i);
        imagex(y_pos(i),   x_pos(i)-1) = imagex(y_pos(i),   x_pos(i)-1) + bit4*energy(i);
        
        % make it proper from all pixel values
        energy(i) = bit0*energy(i)+bit1*energy(i)+bit2*energy(i)+bit3*energy(i)+bit4*energy(i);
        
    end
        
    Cross_Split_Event = table(x_pos,y_pos,energy)
    
end

if (type == 3)
    
    for i = 1:num
        
        choice = randi(3);
        
        if (choice == 1)
       
            imagex = add_split_events(imagex, 1, 2, peak, sigma, 0);

        elseif (choice == 2)
            
            imagex = add_split_events(imagex, 1, 5, peak, sigma, 0);

        else
            
            imagex = add_split_events(imagex, 1, 4, peak, sigma, 0);
            
        end
        
    end
    
end

if (show)
    figure;
    colormap('hot');
    imagesc(imagex);
    colorbar;
end

end