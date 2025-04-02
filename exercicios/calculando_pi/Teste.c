#include <stdio.h>
#include <stdlib.h>

int main(){
    double pi_4 = 0;
    int n = 100;
    int sign = 1;
    int i = 0;
    for (i = 1; i < n; i += 2) {
      if (sign) { pi_4 += 1.0 / i; sign = 0; }
      else { pi_4 -= 1.0 / i; sign = 1; }
    }
    printf("PI = %.12f\n", pi_4 * 4);
}

