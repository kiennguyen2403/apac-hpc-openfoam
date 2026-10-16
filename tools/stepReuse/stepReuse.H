/*---------------------------------------------------------------------------*\
    stepReuse: shared helpers for the LD_PRELOAD replacements in this library
    (gaussSeidelFold.C, kOmegaSSTReuse.C, wallShearStressLazy.C).

    STEP_REUSE_OFF  comma list of features to switch off: gs, komega, wss
                    (the stock code then runs, through the next definition
                    of the replaced symbol)
    STEP_REUSE_REPORT=1  one line per feature on stderr (master) saying
                    which path ran
\*---------------------------------------------------------------------------*/

#ifndef stepReuse_H
#define stepReuse_H

#include <cstdlib>
#include <cstring>
#include <string>

namespace Foam
{
namespace stepReuse
{

inline bool off(const char* key)
{
    const char* e = std::getenv("STEP_REUSE_OFF");
    return e && std::strstr(e, key);
}

//- The next definition (normally the stock one) of a replaced symbol;
//  FatalError when it cannot be found
void* next(const char* mangledName);

//- Print "stepReuse: <msg>" once per key when STEP_REUSE_REPORT is set
void report(const char* key, const std::string& msg);

} // End namespace stepReuse
} // End namespace Foam

#endif
