import numpy as np

def event_consolidator(image, thresh):
    """
    image:  the image to find the peaks in
    thresh: number of sigma to be threshold for events
    
    returns: consolidated, background-subtracted 2D-array
    """
    
    bg    = image.mean()
    sigma = image.std()
    """overall threshold"""
    s = thresh*sigma
    
    image = (image - bg)
    
    """i is the row number"""
    for i in range(len(image)):
        
        """j is the pixel/column number"""
        for j in range(len(image[0])):
            
            pixel = image[i][j]
            
            """if current pixel is bigger than the threshold"""
            if ( pixel > s ):
                image = check2(image,i,j,s,sigma)
    
    return image+bg
    
def check2(image,i,j,s,sigma):
    
    pixel = image[i][j]    
    
    fr = (False)if(i > 0)              else(True)
    lr = (False)if(i < len(image)-1)   else(True)
    fc = (False)if(j > 0)              else(True)
    lc = (False)if(j < len(image[0])-1)else(True)
    
    """surrounding pixels"""
    surr = [[-30000,-30000,-30000],[-30000,-30000,-30000],[-30000,-30000,-30000]]
    
    surr[1][1] = image[i][j]

    """there is a row above it"""
    if (not fr):
        """there is a center above"""
        surr[0][1] = image[i-1][j]
        
        """there is a column before it"""
        if (not fc):
            """there is an upper left"""
            surr[0][0] = image[i-1][j-1]
        
        """there is a column after it"""
        if (not lc):
            """there is an upper right"""
            surr[0][2] = image[i-1][j+1]
    
    """there is a row below it"""
    if (not lr):
        """there is a center below"""
        surr[2][1] = image[i+1][j]
        
        """there is a column before it"""
        if (not fc):
            """there is a lower left"""
            surr[2][0] = image[i+1][j-1]
        
        """there is a column after it"""
        if (not lc):
            """there is a lower right"""
            surr[2][2] = image[i+1][j+1]
    
    """there is a column before it"""
    if (not fc):
        """there is a left"""
        surr[1][0] = image[i][j-1]
    
    """there is a column after it"""
    if (not lc):
        """there is a right"""
        surr[1][2] = image[i][j+1]
    
    surr = np.array(surr)
    
    surr[surr<s] = -30000
    
    """get the maximum value"""
    contact = np.nanmax(surr)
    
    """get its index"""
    y,x = np.where(surr == contact)
    x = x[0]-1
    y = y[0]-1
    
    """
    check all cells
    add the ones that are not BS or the highest
    remove them as they are added
    """
    
    for row in range(len(surr)):
        
        for pix in range(len(surr[row])):
            
            """if it is not BS or the biggest value"""
            if not (surr[row][pix] == -30000) and not (surr[row][pix] == contact):
                
                """corresponding image pixel"""
                ix = pix-1
                iy = row-1
                
                image[i+y][j+x] += image[i+iy][j+ix]
                
                image[i+iy][j+ix] = np.random.normal(0,sigma)
                
    
    return image
    

def check1(image,i,j,s,sigma):
    
    pixel = image[i][j]
    
    """not the first row"""
    if (i > 0):
        up = image[i-1][j]
        
        """the upper pixel"""
        if (up > s):
            """which is larger"""
            if (pixel > up):
                image[i][j]   += up
                image[i-1][j] = np.random.normal(0,sigma)
                
                # image = check(image,i,j,s,sigma)
                
            else:
                image[i-1][j] += pixel
                image[i][j]   = np.random.normal(0,sigma)
                return image
                # image = check(image,i-1,j,s,sigma)
    
    """not the last column"""
    if (i < len(image)-1):
        down = image[i+1][j]
        
        """the lower pixel"""
        if (down > s):
            """which is larger"""
            if (pixel > down):
                image[i][j]   += down
                image[i+1][j] = np.random.normal(0,sigma)
                
                # image = check(image,i,j,s,sigma)
            
            else:
                image[i+1][j] += pixel
                image[i][j]   = np.random.normal(0,sigma)
                return image
                # image = check(image,i+1,j,s,sigma)
    
    """not the first row"""
    if (j > 0):
        left = image[i][j-1]
        
        """the left pixel"""
        if (left > s):
            """which is larger"""
            if (pixel > left):
                image[i][j]   += left
                image[i][j-1] = np.random.normal(0,sigma)
                
                # image = check(image,i,j,s,sigma)
            
            else:
                image[i][j-1] += pixel
                image[i][j]   = np.random.normal(0,sigma)
                return image
                # image = check(image,i,j-1,s,sigma)
    
    """not the last row"""
    if (j < len(image[0])-1):
        right = image[i][j+1]
        
        """the right pixel"""
        if (right > s):
            """which is larger"""
            if (pixel > right):
                image[i][j]   += right
                image[i][j+1] = np.random.normal(0,sigma)
                
                # image = check(image,i,j,s,sigma)
            
            else:
                image[i][j+1] += pixel
                image[i][j]   = np.random.normal(0,sigma)
                return image
                # image = check(image,i,j+1,s,sigma)
    
    if (i > 0):
        if (j > 0):
            
            uleft = image[i-1][j-1]
            
            """upper left pixel"""
            if (uleft > s):
                """which is larger"""
                if (pixel > uleft):
                    image[i][j]     += uleft
                    image[i-1][j-1] = np.random.normal(0,sigma)
                    
                    # image = check(image,i,j,s,sigma)
                    
                else:
                    image[i-1][j-1] += pixel
                    image[i][j]     = np.random.normal(0,sigma)
                    return image
                    # image = check(image,i-1,j-1,s,sigma)
        
        if (j < len(image[0])-1):
            
            uright = image[i-1][j+1]
            
            """upper right"""
            if (uright > s):
                """which is larger"""
                if (pixel > uright):
                    image[i][j]     += uright
                    image[i-1][j+1] = np.random.normal(0,sigma)
                    
                    # image = check(image,i,j,s,sigma)
                
                else:
                    image[i-1][j+1] += pixel
                    image[i][j]     = np.random.normal(0,sigma)
                    return image
                    # image = check(image,i-1,j+1,s,sigma)
    
    if (i < len(image)-1):
        if (j > 0):
            
            dleft = image[i+1][j-1]
            
            """lower left"""
            if (dleft > s):
                """which is larger"""
                if (pixel > dleft):
                    image[i][j]     += dleft
                    image[i+1][j-1] = np.random.normal(0,sigma)
                    
                    # image = check(image,i,j,s,sigma)
                
                else:
                    image[i+1][j-1] += pixel
                    image[i][j]     = np.random.normal(0,sigma)
                    return image
                    # image = check(image,i+1,j-1,s,sigma)
        
        if (j < len(image[0])-1):
            
            dright = image[i+1][j+1]
            
            """lower right"""
            if (dright > s):
                """which is larger"""
                if (pixel > dright):
                    image[i][j]     += dright
                    image[i+1][j+1] = np.random.normal(0,sigma)
                    
                    # image = check(image,i,j,s,sigma)
                
                else:
                    image[i+1][j+1] += pixel
                    image[i][j]     = np.random.normal(0,sigma)
                    return image
                    # image = check(image,i+1,j+1,s,sigma)
    
    return image
    
    
    