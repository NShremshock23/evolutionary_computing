#include <stdio.h>
#include <stdlib.h>

int main() {
    char a = '0';
    char b = '1';

    char c = '0' + ( (int)a ^ (int)b );


    printf("%c", c);
    return 0;
}
