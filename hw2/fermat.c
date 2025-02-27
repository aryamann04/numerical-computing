//
//  fermat.c
//
//
//  Created by Aryaman Nagpal on 2/2/25.
//

#include <stdio.h>
#include <math.h>

void fermat_float(int n, float x, float y, float z);
void fermat_double(int n, double x, double y, double z);

int main() {
    
    // part (a)
    fermat_float(4, 717, 967, 1033);
    
    // part (b)
    fermat_double(4, 717, 967, 1033);
    
    // part (c)
    fermat_float(5, 844487, 1288439, 1318202);
    fermat_double(5, 844487, 1288439, 1318202);
    
    return 0;
}

void fermat_float(int n, float x, float y, float z) {
    
    float a = pow(x, n) + pow(y, n);
    float b = pow(z, n);
    
    printf("x^n + y^n = %13.7e\n", a);
    printf("z^n = %13.7e\n", b);
    printf("%d\n\n", a == b);
}

void fermat_double(int n, double x, double y, double z) {
    
    double a = pow(x, n) + pow(y, n);
    double b = pow(z, n);
    
    printf("x^n + y^n = %22.16e\n", a);
    printf("z^n = %22.16e\n", b);
    printf("%d\n\n", a == b);
}
