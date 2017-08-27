
#ifndef BL_BARENA_H
#define BL_BARENA_H

#include <AMReX_Arena.H>


namespace amrex {
/**
* \brief A Concrete Class for Dynamic Memory Management
* This is the simplest dynamic memory management class derived from Arena.
* Makes calls to ::operator new() and ::operator delete().
*
* This class can allocate 3 types of dynamic memory:
* 1) CPU pageable memory: alloc(...), free(...)
* 2) CPU unpageable (pinned) memory: alloc_pinned(...), free_pinned(...)
* 3) GPU memory: alloc_device(...), free_device(...)
*/

class BArena
    :
    public Arena
{
public:
    /**
    * \brief Allocates a dynamic memory arena of size sz.
    * Returns a pointer to this memory.
    * This implementation allocate pageable CPU memory
    */
    virtual void* alloc (std::size_t sz) override;

    /**
    * \brief Allocates a dynamic memory arena of size sz.
    * Returns a pointer to this memory.
    * This implementation allocate unpageable CPU memory
    * It has the same effect as alloc if CUDA is not defined
    */
    virtual void* alloc_pinned (std::size_t sz) override;

    /**
    * \brief Deletes the arena pointed to by pt.
    */
    virtual void free (void* pt) override;

    /**
    * \brief Deletes the arena pointed to by pt.
    * Used to free pinned memory
    */
    virtual void free_pinned (void* pt) override;

    /**
    * \brief Allocates a dynamic device memory arena of size sz.
    * Returns a pointer to this memory.
    */
    virtual void* alloc_device (std::size_t sz, int device_id) override;
#ifdef CUDA_ARRAY
    // for 2d array on GPU
    virtual void* alloc_device_2d (std::size_t& _pitch, std::size_t _isize, std::size_t _jsize, int device_id) override;
#endif
    /**
    * \brief Deletes the device arena pointed to by pt.
    */
    virtual void free_device (void* pt, int device_id) override;
};

}

#endif /*BL_BARENA_H*/
