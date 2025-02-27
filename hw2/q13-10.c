//
//  q13-10.c
//  
//
//  Created by Aryaman Nagpal on 2/8/25.
//

#include <stdio.h>
#include <math.h>

void program_13_3(double x);

int main() {
    
    double x = expm1(0.000001);
    printf("part (ii): %22.16e\n\n", x);
    
    double y = exp(0.000001) - 1;
    printf("part (iii): %22.16e\n\n", y);
    puts("part (iv):");
    program_13_3(0.000001);
}

void program_13_3(double x) {
    int n = 1;  // start at n = 1 instead
    double term = x, oldsum = 0.0, newsum = x; // begin from +x instead of +1
    
    /* terminates when the new sum is no different from the old sum */
    while (newsum != oldsum){
        oldsum = newsum;
        n++;
        term = (term*x)/n; /* term has the value (x^n)/(n!) */
        newsum = newsum + term; /* approximates exp(x) */
        printf("n = %3d term = %13.6e newsum = %13.6e \n",
        n,term,newsum);
    }
    
    printf("From summing the series, expm1(x)=%22.16e \n", newsum);
    printf("Using the standard function, expm1(x)=%22.16e \n", expm1(x));
}
