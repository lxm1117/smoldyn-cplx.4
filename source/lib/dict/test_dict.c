#include <stdio.h>
#include <assert.h>

#include "dict.h"

int
main()
{
    Dict d;
    char buf[512];
    int i;

    d = DictCreate();

    DictInsert(d, "foo", "hello");		// h= (multiplier*"f"+"o")*multiplier+"o" 
    DictInsert(d, "foo", "hello2");
	DictInsert(d, "koo", "test");
	DictInsert(d, "foo", "hello");
	DictInsert(d, "foo", "hello3");
	DictInsert(d, "koo", "test1");
	DictInsert(d, "gii", "latest");
    DictDelete(d, "foo");
    DictDelete(d, "foo");
    assert(DictSearch(d, "foo") == 0);
    DictDelete(d, "foo");

    for(i = 0; i < 10000; i++) {
        sprintf(buf, "%d", i);
        DictInsert(d, buf, buf);
    }

    DictDestroy(d);
    return 0;
}

    
