#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "dict1.h"

struct elt {
    struct elt *next;
    int key1;
    int key2;
	int value;	
};

struct dict {
    int size;           /* size of the pointer table */
    int n;              /* number of elements stored */
    struct elt **table;
};

#define INITIAL_SIZE (2046)
#define GROWTH_FACTOR (2)
#define MAX_LOAD_FACTOR (1)

/* dictionary initialization code used in both DictCreate and grow */
Dict internalDictCreate(int size)
{
    Dict d;
    int i;

    d = malloc(sizeof(*d));

    assert(d != 0);

    d->size = size;
    d->n = 0;
    d->table = malloc(sizeof(struct elt *) * d->size);

    assert(d->table != 0);

    for(i = 0; i < d->size; i++) d->table[i] = 0;

    return d;
}

Dict DictCreate(void)
{
    return internalDictCreate(INITIAL_SIZE);
}

void DictDestroy(Dict d)
{
    int i;
    struct elt *e;
    struct elt *next;

    for(i = 0; i < d->size; i++) {
        for(e = d->table[i]; e != 0; e = next) {
            next = e->next;

            // free(e->key1);
			// free(e->key2);
            // free(e->value);
            free(e);
        }
    }

    free(d->table);
    free(d);
}

#define MULTIPLIER (97)
/* stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
	Szudzik's paring	*/

static unsigned long hash_function(const int a, const int b)
{
	unsigned long A;
	unsigned long B;
	unsigned long h;

	A= a>0? 2*(long)a:-2*(long)a-1;
	B= b>0? 2*(long)b:-2*(long)b-1;	
	h= (A>=B? A*A+A+B:A+B*B)/2; 	
	h= a>=0 && b>=0? h: -h-1;
    return h;
}

static void grow(Dict d)
{
    Dict d2;            /* new dictionary we'll create */
    struct dict swap;   /* temporary structure for brain transplant */
    int i;
    struct elt *e;

    d2 = internalDictCreate(d->size * GROWTH_FACTOR);

    for(i = 0; i < d->size; i++) {
        for(e = d->table[i]; e != 0; e = e->next) {
            /* note: this recopies everything */
            /* a more efficient implementation would
             * patch out the strdups inside DictInsert
             * to avoid this problem */
            DictInsert(d2, e->key1, e->key2, e->value);
        }
    }

    /* the hideous part */
    /* We'll swap the guts of d and d2 */
    /* then call DictDestroy on d2 */
    swap = *d;
    *d = *d2;
    *d2 = swap;

    DictDestroy(d2);
}

/* insert a new key-value pair into an existing dictionary */
void DictInsert(Dict d, const int key1, const int key2, const int value)
{
    struct elt *e;
    unsigned long h;

    assert(key1);
	assert(key2);
    assert(value);

    e = malloc(sizeof(*e));

    assert(e);

    e->key1 = key1;
	e->key2 = key2;
    e->value = value;

    h = hash_function(key1,key2) % d->size;

    e->next = d->table[h];
    d->table[h] = e;

    d->n++;

    /* grow table if there is not enough room */
    if(d->n >= d->size * MAX_LOAD_FACTOR) {
        grow(d);
    }
}

/* return the most recently inserted value associated with a key */
/* or 0 if no matching key is present */
const int DictSearch(Dict d, const int key1, const int key2)
{
    struct elt *e;

    for(e = d->table[hash_function(key1, key2) % d->size]; e != 0; e = e->next) {
        // if(!strcmp(e->key, key)) {
            /* got it */
        //    return e->value;
        // }
		if(e->key1==key1 && e->key2==key2) return e->value;
	
    }

    return 0;
}


void DictDelete(Dict d, const int key1, const int key2)
{
	struct elt **prev;
	struct elt *e;

	for(prev= &(d->table[hash_function(key1,key2) % d->size]); *prev!=0; prev=&((*prev)->next)){
		if((*prev)->key1=key1, (*prev)->key2=key2){
			e = *prev;	
			*prev = e->next;
			free(e);
			return;
		}	
	}
}
