
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>

namespace amrex
{

struct FabInfo
{
    bool alloc = true;
    bool device = false;
    bool shared = false;
#ifdef CUDA
    int device_id = 0;
#endif

    FabInfo& SetAlloc (bool a) {
        alloc = a;
        return *this;
    }

    FabInfo& SetDevice (bool d) {
	device = d;
	return *this;
    }

    FabInfo& SetShared (bool s) {
        shared = s;
        return *this;
    }

#ifdef CUDA
    FabInfo& useDevice (int id) {
	device_id = id;
	return *this;
    }
#endif
};

template <class FAB>
class FabFactory
{
public:
    virtual ~FabFactory () {}
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const = 0;
};

template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
    {
        FAB* fab = new FAB(box, ncomps, info.alloc, info.shared);
#ifdef CUDA
	if (info.device) fab->allocate_device(info.device_id);
#endif
	return fab;
    }
};

}

#endif
