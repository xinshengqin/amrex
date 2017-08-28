#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <limits>
#include <algorithm>
#include <new>
#include <memory>
#include <cstring>
#include <cstdint>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_CArena.H>
#ifdef CUDA
#include <AMReX_DArena.H>
#include <AMReX_EArena.H>
#include <AMReX_CUDA_helper.H>
#include <cuda_runtime_api.h>
#endif
#include <AMReX_MemPool.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifndef FORTRAN_BOXLIB
#include <AMReX_ParmParse.H>
#endif

using namespace amrex;

namespace
{
    static Array<std::unique_ptr<CArena> > the_memory_pool;
#ifdef CUDA
    // Each DArena manages memory in one GPU
    // device_memory_pool[i] manages memory in GPU i
    static Array<std::unique_ptr<DArena> > device_memory_pool;

    // EArena manages pinned (unpageable) CPU memory for asynchronous 
    // data transfer between CPU and GPU
    static std::unique_ptr<EArena> pinned_memory_pool;
#endif
#if defined(BL_TESTING) || defined(DEBUG)
    static int init_snan = 1;
#else
    static int init_snan = 0;
#endif
}

extern "C" {

void amrex_mempool_init()
{
    static bool initialized = false;
    if (!initialized)
    {
	initialized = true;

#ifndef FORTRAN_BOXLIB
        ParmParse pp("fab");
	pp.query("init_snan", init_snan);
#endif

#ifdef _OPENMP
	int nthreads = omp_get_max_threads();
#else
	int nthreads = 1;
#endif
	the_memory_pool.resize(nthreads);
	for (int i=0; i<nthreads; ++i) {
	    the_memory_pool[i].reset(new CArena);
	}
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    size_t N = 1024*1024*sizeof(double);
	    void *p = amrex_mempool_alloc(N);
	    memset(p, 0, N);
	    amrex_mempool_free(p);
	}

#ifdef BL_MEM_PROFILING
	MemProfiler::add("MemPool", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     int MB_min, MB_max, MB_tot;
			     amrex_mempool_get_stats(MB_min, MB_max, MB_tot);
			     long b = MB_tot * (1024L*1024L);
			     return {b, b};
			 }));
#endif
#ifdef CUDA
        int num_devices = ParallelDescriptor::get_num_devices_used();
	device_memory_pool.resize(num_devices);
	for (int i=0; i<num_devices; ++i) {
	    device_memory_pool[i].reset(new DArena(i,0));
	}
        pinned_memory_pool.reset(new EArena);
#endif
    }
}

void* amrex_mempool_alloc (size_t nbytes)
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return the_memory_pool[tid]->alloc(nbytes);
}


void amrex_mempool_free (void* p) 
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  the_memory_pool[tid]->free(p);
}

#ifdef CUDA
void* amrex_mempool_alloc_pinned (size_t nbytes)
{
  return pinned_memory_pool->alloc_pinned(nbytes);
}

void amrex_mempool_free_pinned (void* p) 
{
  pinned_memory_pool->free_pinned(p);
}

void* amrex_mempool_alloc_gpu (size_t nbytes, int device_id)
{

    BL_ASSERT(device_id >= 0);
    return device_memory_pool[device_id]->alloc_device(nbytes, device_id);
}

void amrex_mempool_free_gpu (void* p, int device_id) 
{
    BL_ASSERT(device_id >= 0);
    device_memory_pool[device_id]->free_device(p, device_id);
}
#endif

void amrex_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot) // min, max & tot in MB
{
  size_t hsu_min=std::numeric_limits<size_t>::max();
  size_t hsu_max=0;
  size_t hsu_tot=0;
  for (const auto& mp : the_memory_pool) {
    size_t hsu = mp->heap_space_used();
    hsu_min = std::min(hsu, hsu_min);
    hsu_max = std::max(hsu, hsu_max);
    hsu_tot += hsu;
  }
  mp_min = hsu_min/(1024*1024);
  mp_max = hsu_max/(1024*1024);
  mp_tot = hsu_tot/(1024*1024);
}

// We should eventually use Real instead of double.
// We cannot do it now because of F_BaseLib.
void amrex_real_array_init (double* p, size_t nelems)
{
    if (init_snan) amrex_array_init_snan(p, nelems);
}

void amrex_array_init_snan (double* p, size_t nelems)
{
#ifdef UINT64_MAX
    const uint64_t snan = UINT64_C(0x7ff0000080000001);
#else
    static_assert(sizeof(double) == sizeof(long long), "MemPool: sizeof double != sizeof long long");
    const long long snan = 0x7ff0000080000001LL;
#endif
    for (size_t i = 0; i < nelems; ++i) {
        std::memcpy(p++, &snan, sizeof(double));
    }
}

}
