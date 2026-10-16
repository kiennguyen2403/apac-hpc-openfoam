#include "stepReuse.H"
#include "error.H"
#include "UPstream.H"

#include <dlfcn.h>
#include <iostream>
#include <set>

void* Foam::stepReuse::next(const char* mangledName)
{
    void* p = ::dlsym(RTLD_NEXT, mangledName);

    if (!p)
    {
        FatalErrorInFunction
            << "stepReuse: no stock definition of " << mangledName
            << exit(FatalError);
    }

    return p;
}


void Foam::stepReuse::report(const char* key, const std::string& msg)
{
    static std::set<std::string> done;

    if (!std::getenv("STEP_REUSE_REPORT") || !UPstream::master())
    {
        return;
    }

    if (done.insert(key).second)
    {
        std::cerr << "stepReuse: " << msg << std::endl;
    }
}
