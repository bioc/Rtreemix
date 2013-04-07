#include <Rdefines.h>
#include "include/Rtreemix_patch.h"

void _Rtreemix_exit(int status)
{
    Rf_error
        ("internal: mtreemix invoked 'exit(%d)'; see warnings() and restart R",
         status);
}
