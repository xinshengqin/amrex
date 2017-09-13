#ifndef _TINY_PROFILER_H_
#define _TINY_PROFILER_H_

#include <string>
#include <stack>
#include <map>
#include <utility>
#include <limits>

#include <AMReX_REAL.H>

#ifdef CUDA
#include "nvToolsExt.h"
#endif


namespace amrex {

#ifdef CUDA
namespace nvtx {
extern uint32_t colors[];
extern int num_colors;
#define PUSH_RANGE(name,cid) { \
        int color_id = cid; \
        color_id = color_id%num_colors;\
        nvtxEventAttributes_t eventAttrib = {0}; \
        eventAttrib.version = NVTX_VERSION; \
        eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
        eventAttrib.colorType = NVTX_COLOR_ARGB; \
        eventAttrib.color = colors[color_id]; \
        eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
        eventAttrib.message.ascii = name; \
        nvtxRangePushEx(&eventAttrib); \
}
#define POP_RANGE nvtxRangePop();

class CudaTracer {
    public:
        CudaTracer(const char* name, int cid = 0) {
            PUSH_RANGE(name,cid);
        }
        ~CudaTracer() {
            POP_RANGE;
        }
};

}
#endif


//! A simple profiler that returns basic performance information (e.g. min, max, and average running time)
class TinyProfiler
{
public:
    TinyProfiler (const std::string &funcname);
    TinyProfiler (const std::string &funcname, bool start_);
    ~TinyProfiler ();

    void start ();
    void stop ();

    static void Initialize ();
    static void Finalize ();

private:
    //! stats on a single process
    struct Stats   
    {
	Stats () : depth(0), n(0L), dtin(0.0), dtex(0.0) { } 
	int  depth; // recursive depth
	long n;     // number of calls
	Real dtin;  // inclusive dt
	Real dtex;  // exclusive dt
    };

    //! stats across processes
    struct ProcStats 
    {
	ProcStats () : nmin(std::numeric_limits<long>::max()), 
		       navg(0L), nmax(0L), 
		       dtinmin(std::numeric_limits<Real>::max()), 
		       dtinavg(0.0), dtinmax(0.0),
		       dtexmin(std::numeric_limits<Real>::max()), 
		       dtexavg(0.0), dtexmax(0.0)  {}
	long nmin, navg, nmax;
	Real dtinmin, dtinavg, dtinmax;
	Real dtexmin, dtexavg, dtexmax;
	std::string fname;
	static bool compex (const ProcStats& lhs, const ProcStats& rhs) {
	    return lhs.dtexmax > rhs.dtexmax;
	}
	static bool compin (const ProcStats& lhs, const ProcStats& rhs) {
	    return lhs.dtinmax > rhs.dtinmax;
	}
    };

    std::string fname;
    bool running;
    int global_depth;

    static std::stack<std::pair<Real,Real> > ttstack;
    static std::map<std::string, Stats> statsmap;
    static Real t_init;

#ifdef CUDA
    nvtxRangeId_t nvtx_id;
#endif
};

}

#endif
