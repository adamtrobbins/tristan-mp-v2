#define STR_MAX   280
#define END_LINE  char(10)
#define TAB       char(9)

#define M_PI      3.141592653589793
#define CC        0.45

#ifdef threeD
#define fldBoundZ           (-NGHOST) : ((this_meshblock%ptr%sz) - 1 + (NGHOST))
#define globalZ             (global_mesh%sz)
#define globalZwithGhosts   ((global_mesh%sz) + ((2) * (NGHOST)))
#else
#define fldBoundZ           (0) : (0)
#define globalZ             (1)
#define globalZwithGhosts   (1)
#endif
