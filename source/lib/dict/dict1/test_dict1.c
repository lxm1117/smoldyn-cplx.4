#include <stdio.h>
// #include <assert.h>

#include "dict1.h"

int main()
{
    Dict d;
    // char buf[512];
    d = DictCreate();

    DictInsert(d, 1, 10, 1);		// "ca", "cam", 1
	DictInsert(d, 1, 10, 5);
	DictInsert(d, 1, 10, 21);
	DictInsert(d, 10, 1, 4);		// "cam", "ca", 4
	DictInsert(d, 1,  12, 1);		// "cam", "bK", 1 
	DictInsert(d, 12, 1 ,2);		// "bK", "cam", 2
	
	DictSearch(d, 1, 10);
	DictSearch(d, 10, 1);

    DictDelete(d, 12, 1);
    DictSearch(d, 12, 1);
    // assert(DictSearch(d, "foo") == 0);

    // for(i = 0; i < 10000; i++) {
    //    sprintf(buf, "%d", i);
    //    DictInsert(d, buf, buf);
    // }

    DictDestroy(d);
    return 0;
}

    
