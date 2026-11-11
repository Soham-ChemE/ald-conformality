/*---------------------------------------------------------------------------*\
    aldkinHandle: owning wrapper of an aldkin_t handle (lib/aldkin) for the OpenFOAM adapters.

    Copying clones the kinetics state (coverages, faces, profiles); the destructor releases it. Every library
    return code is checked by check(), which raises a FatalError carrying aldkin_last_error().
\*---------------------------------------------------------------------------*/
#ifndef aldkinHandle_H
#define aldkinHandle_H

#include "aldkin.h"
#include "error.H"
#include "scalarField.H"

namespace Foam
{

class aldkinHandle
{
    aldkin_t* h_;

public:

    aldkinHandle() : h_(nullptr) {}

    explicit aldkinHandle(const aldkin_params& p)
    :
        h_(aldkin_create(&p))
    {
        if (!h_)
        {
            FatalErrorInFunction << "aldkin_create failed: " << aldkin_last_error(nullptr) << exit(FatalError);
        }
    }

    aldkinHandle(const aldkinHandle& other)
    :
        h_(other.h_ ? aldkin_clone(other.h_) : nullptr)
    {
        if (other.h_ && !h_)
        {
            FatalErrorInFunction << "aldkin_clone failed: " << aldkin_last_error(nullptr) << exit(FatalError);
        }
    }

    aldkinHandle& operator=(const aldkinHandle& other)
    {
        if (this != &other)
        {
            aldkin_t* n = other.h_ ? aldkin_clone(other.h_) : nullptr;
            if (other.h_ && !n)
            {
                FatalErrorInFunction << "aldkin_clone failed: " << aldkin_last_error(nullptr) << exit(FatalError);
            }
            aldkin_destroy(h_);
            h_ = n;
        }
        return *this;
    }

    ~aldkinHandle() { aldkin_destroy(h_); }

    aldkin_t* get() const { return h_; }
    bool valid() const { return h_ != nullptr; }

    void check(int rc, const char* what) const
    {
        if (rc != ALDKIN_OK)
        {
            FatalErrorInFunction << what << " failed (code " << rc << "): " << aldkin_last_error(h_) << exit(FatalError);
        }
    }

    // coverage field of one selector (n values)
    scalarField theta(int which, label n) const
    {
        scalarField f(n, 0.0);
        check(aldkin_theta(h_, which, f.data()), "aldkin_theta");
        return f;
    }
};

} // End namespace Foam

#endif
