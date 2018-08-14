function [ image ] = surrounding( imagex, i, j, s, sigma )

[ Y, X ] = size(imagex);

pixel = imagex(i,j);

if(i>1);fr=false;else;fr=true;end;
if(i<Y);lr=false;else;lr=true;end;
if(j>1);fc=false;else;fc=true;end;
if(j<X);lc=false;else;lc=true;end;

surr =[[-30000,-30000,-30000];[-30000,-30000,-30000];[-30000,-30000,-30000]];

surr(2,2) = imagex(i,j);

% there is a row above it
if (~ fr)
    % there is a center above
    surr(1,2) = imagex(i-1,j);

    % there is a column before it
    if (~ fc)
        % there is an upper left
        surr(1,1) = imagex(i-1,j-1);
    end

    % there is a column after it
    if (~ lc)
        % there is an upper right
        surr(1,3) = imagex(i-1,j+1);
    end
end
% there is a row below it
if (~ lr)
    % there is a center below
    surr(3,2) = imagex(i+1,j);

    % there is a column before it
    if (~ fc)
        % there is a lower left
        surr(3,1) = imagex(i+1,j-1);
    end

    % there is a column after it
    if (~ lc)
        % there is a lower right
        surr(3,3) = imagex(i+1,j+1);
    end
end
% there is a column before it
if (~ fc)
    % there is a left
    surr(2,1) = imagex(i,j-1);
end
% there is a column after it
if (~ lc)
    % there is a right
    surr(2,3) = imagex(i,j+1);
end

image = imagex;
surr(surr < s) = -30000;
contact = max(max(surr));
[ y, x ] = find(surr == contact);
x = x-2;
y = y-2;

for row = 1:3
    
    for pix = 1:3
        
        if (~(surr(row,pix)==-30000))&(~(surr(row,pix)==contact))
            
            ix = pix-2;
            iy = row-2;
            
            image(i+y,j+x) = image(i+y,j+x) + image(i+iy,j+ix);
            
            image(i+iy,j+ix) = random('norm',0,sigma);
            
        end
        
    end
    
end

end