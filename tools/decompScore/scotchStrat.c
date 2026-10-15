/* Print a scotch 6.1.0 graph-mapping strategy built from flags, and check it parses back.
   usage: scotchStrat <flags> <nParts> <imbalance> */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "scotch.h"
int main(int argc, char **argv)
{
    if (argc != 4) { fprintf(stderr, "usage: %s flags nParts imbalance\n", argv[0]); return 2; }
    SCOTCH_Strat s, t;
    char buf[1 << 16];
    FILE *m = fmemopen(buf, sizeof buf, "w");
    SCOTCH_stratInit(&s);
    if (SCOTCH_stratGraphMapBuild(&s, (SCOTCH_Num) strtol(argv[1], 0, 0), (SCOTCH_Num) atoi(argv[2]), atof(argv[3]))) return 1;
    SCOTCH_stratSave(&s, m);
    fclose(m);
    SCOTCH_stratInit(&t);
    if (SCOTCH_stratGraphMap(&t, buf)) { fprintf(stderr, "reparse failed\n"); return 1; }
    puts(buf);
    return 0;
}
