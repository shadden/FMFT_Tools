#include <stdio.h>

int main(int argc, char *argv[]){
 int i; int tmp;
 for(i=0; i < argc; i++){
   printf("argv[%d] is: %s\n",i,argv[i]);
 }
 sscanf(argv[2],"%d",&tmp);
 printf(" %d + 10 = %d\n",tmp,tmp+10);
 return 1;
}
