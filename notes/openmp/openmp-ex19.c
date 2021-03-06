#include <stdio.h>
#include <unistd.h>
#include <omp.h>

int main(void)
{
  #pragma omp parallel
  {
    int thread_num = omp_get_thread_num(), i;
    char wait[BUFSIZ] = {'\0'};
    for (i = 0; i < 4 * thread_num; i++) wait[i] = ' ';
    sleep(thread_num);
    printf("%srow row row your boat...\n",wait);
    /* But we can add synchronization points all over the place */
    #pragma omp barrier
    sleep(thread_num);
    printf("%s...gently down the stream...\n",wait);
  }
  printf("Better.\n");
  return 0;
}
