//gcc -o main2 Main2.c -lpthread

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define N 1000
#define THREADS 2

typedef struct node {
    int value;
    struct node *next;
} node_t;

typedef struct {
    node_t *head;
    pthread_mutex_t lock;
} list_t;

list_t *lists;
int K;

void insert(list_t *list, int value) {
    pthread_mutex_lock(&list->lock);
    node_t *new_node = malloc(sizeof(node_t));
    new_node->value = value;
    new_node->next = list->head;
    list->head = new_node;
    pthread_mutex_unlock(&list->lock);
}

void *thread_func(void *arg) {
    int i;
    for (i = 0; i < N; i++) {
        int value = rand();
        int which = rand() % K;
        insert(&lists[which], value);
    }
    return NULL;
}

int main() {
    printf("Insira o número de listas: ");
    scanf("%d", &K);

    lists = malloc(K * sizeof(list_t));
    for (int i = 0; i < K; i++) {
        lists[i].head = NULL;
        pthread_mutex_init(&lists[i].lock, NULL);
    }

    srand(time(NULL));
    pthread_t threads[THREADS];
    for (int i = 0; i < THREADS; i++)
        pthread_create(&threads[i], NULL, thread_func, NULL);
    for (int i = 0; i < THREADS; i++)
        pthread_join(threads[i], NULL);

    for (int i = 0; i < K; i++) {
        int count = 0;
        node_t *cur;
        for (cur = lists[i].head; cur; cur = cur->next) count++;
        printf("Lista %d: %d nós\n", i, count);
    }

    return 0;
}