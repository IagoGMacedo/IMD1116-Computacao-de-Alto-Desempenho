//gcc -o main1 Main1.c -lpthread

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define N 1000

typedef struct node {
    int value;
    struct node *next;
} node_t;

typedef struct {
    node_t *head;
    pthread_mutex_t lock;
} list_t;

list_t list1 = {NULL, PTHREAD_MUTEX_INITIALIZER};
list_t list2 = {NULL, PTHREAD_MUTEX_INITIALIZER};

void insert(list_t *list, int value) {
    pthread_mutex_lock(&list->lock); // Região crítica nomeada
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
        int which = rand() % 2;
        if (which == 0)
            insert(&list1, value);
        else
            insert(&list2, value);
    }
    return NULL;
}

int main() {
    srand(time(NULL));
    pthread_t t1, t2;
    pthread_create(&t1, NULL, thread_func, NULL);
    pthread_create(&t2, NULL, thread_func, NULL);
    pthread_join(t1, NULL);
    pthread_join(t2, NULL);

    int count1 = 0, count2 = 0;
    node_t *cur;
    for (cur = list1.head; cur; cur = cur->next) count1++;
    for (cur = list2.head; cur; cur = cur->next) count2++;
    printf("List1: %d nodes\nList2: %d nodes\n", count1, count2);
    return 0;
}