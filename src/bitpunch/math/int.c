/*
This file is part of BitPunch
Copyright (C) 2013-2015 Frantisek Uhrecky <frantisek.uhrecky[what here]gmail.com>
Copyright (C) 2013-2014 Andrej Gulyas <andrej.guly[what here]gmail.com>
Copyright (C) 2013-2014 Marek Klein  <kleinmrk[what here]gmail.com>
Copyright (C) 2013-2014 Filip Machovec  <filipmachovec[what here]yahoo.com>
Copyright (C) 2013-2014 Jozef Kudlac <jozef[what here]kudlac.sk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <math.h>
#include "int.h"

int BPU_isPrime(int n) {
    int i;

    if (n == 2)
        return 1;
    for (i = 2; i <= sqrt((double) n); i++) {
        if (n % i == 0)
            return 0;
    }
    return 1;
}

int msb32(uint32_t x, int start, int len, int ele_size) {
    // crop lenght to half
    len /= 2;

    // last recursion
    if (len == 1) {
        // return position of highest bit
        if ((x >> len) == 1ul)
            return start + 1;
        else
            return start;
    }

    // right half of integer
    if ((x >> len) != 0ul) {
        start += len;
        x >>= len;
    }
    // left half of integer
    else {
        x <<= ele_size - len;
        x >>= ele_size - len;
    }

    // recursion
    return msb32(x, start, len, ele_size);
}

int min(int x, int y) {
    
    if(x <= y){
        return x;
    }
    
    return y;
}

int absVal(int x) {
    
    if(x < 0){
        return -x;
    }
    
    return x;
}

/*long long int binomialCoeff(int n, int k){ 
  // Base Cases 
  if (k == 0 || k == n){
    return 1; 
  }
  
  // Recur 
  return  (binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k)); 
} */


#define MAX 100
 
long long int binomialCoeff(int m, int n)
{
    long long int mat[MAX][MAX];
    int i, j;
    if (n > m) return 0;
    if( (n == 0) || (m == n) ) return 1;
    for(j = 0; j < n; j ++)
    {
        mat[0][j] = 1;
        if (j == 0)
        {
            for (i = 1; i<= m - n; i ++ ) mat[i][j] = i + 1;
        }
        else
        {
            for (i = 1; i<= m - n; i ++) mat[i][j] = mat[i - 1][j] + mat[i][j - 1];
        }
    }
    return (mat[m - n][n - 1]);
}