//
//  central-difference.c
//  
//
//  Created by Aryaman Nagpal on 2/2/25.
//

#include <stdio.h>
#include <math.h>

int main() {
    
    int n = 1;
    double x = 1.0, h = 1.0, deriv = cos(x), centraldiff, error;
    double minh, minerror = 1;
    
    printf(" deriv = %13.6e \n", deriv);
    printf("h        centraldiff   abs(deriv - centraldiff) \n");
    
    while(n <= 20) {
        
        h = h / 10;
        centraldiff = (sin(x+h) - sin(x-h))/(2 * h);
        error = fabs(deriv - centraldiff);
        printf("%5.1e %13.6e %13.6e \n", h, centraldiff, error);
        
        if (error < minerror) {
            minh = h;
            minerror = error;
        }
        
        n++;
    }
    
    printf("min error is %13.6e achieved at h = %5.1e", minerror, minh);
    
    return 0;
}


