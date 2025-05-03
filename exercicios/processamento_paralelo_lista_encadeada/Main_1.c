// gcc -fopenmp -o Main Main.c

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

typedef struct Node
{
    char filename[100];
    struct Node *next;
} Node;

Node *create_node(const char *filename)
{
    Node *new_node = (Node *)malloc(sizeof(Node));
    snprintf(new_node->filename, 100, "%s", filename);
    new_node->next = NULL;
    return new_node;
}

void free_list(Node *head)
{
    Node *tmp;
    while (head)
    {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

int main()
{
    Node *head = create_node("file1.txt");
    head->next = create_node("file2.txt");
    head->next->next = create_node("file3.txt");
    head->next->next->next = create_node("file4.txt");

    Node *curr = head;

#pragma omp parallel
    {

        while (curr)
        {
            Node *node = curr;
#pragma omp task
            {
                printf("Arquivo: %s | Thread: %d\n", node->filename, omp_get_thread_num());
            }
            curr = curr->next;
        }
    }

    free_list(head);
    return 0;
}