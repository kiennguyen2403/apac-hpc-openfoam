#include "fixedIterSolvers.H"
#include "error.H"
#include "UPstream.H"

#include <dlfcn.h>
#include <iostream>
#include <set>

void* Foam::fixedIterSolvers::next(const char* mangledName)
{
    void* p = ::dlsym(RTLD_NEXT, mangledName);

    if (!p)
    {
        FatalErrorInFunction
            << "fixedIterSolvers: no stock definition of " << mangledName
            << exit(FatalError);
    }

    return p;
}


bool Foam::fixedIterSolvers::reporting()
{
    static const bool on = std::getenv("STEP_REUSE_REPORT") != nullptr;

    return on && UPstream::master();
}


void Foam::fixedIterSolvers::report(const std::string& msg)
{
    static std::set<std::string> done;

    if (reporting() && done.insert(msg).second)
    {
        std::cerr << "fixedIterSolvers: " << msg << std::endl;
    }
}
