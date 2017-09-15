#include <iostream>
#include <array>
#include <vector>
#include <utility>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <AmrCoreAdv.H>
#include <AmrCoreAdvPhysBC.H>
#include <AmrCoreAdv_F.H>

#ifdef CUDA
#include <cuda_runtime_api.h>
#include <AMReX_MemPool.H>
#include <cublas_v2.h>
#endif

using namespace amrex;

static const bool my_verbose = false;

#ifdef CUDA
struct CopyTag {
    public:
        std::array<Box, BL_SPACEDIM> box;
        std::array<FArrayBox*, BL_SPACEDIM> dst;
        std::array<FArrayBox*, BL_SPACEDIM> src;
        CopyTag () {}
        // copy constructor
        CopyTag (const CopyTag& rhs) {
            std::cout << "copy constructor called." << std::endl;
            for (int i = 0; i < BL_SPACEDIM ; i++) {
                box[i] = rhs.box[i];
                dst[i] = rhs.dst[i];
                src[i] = rhs.src[i];
            }
        }
};

void CUDART_CB cudaCallback_copyFAB(cudaStream_t stream, cudaError_t status, void *data){
    checkCudaErrors(status);
    // CopyTag* copy_tag_pt = (CopyTag*) data;
    CopyTag& copy_tag = *((CopyTag*) data);
    for (int i = 0; i < BL_SPACEDIM ; i++) {
        (copy_tag.dst[i])->copy(*(copy_tag.src[i]), copy_tag.box[i]);
    }
}
#endif

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
AmrCoreAdv::AmrCoreAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually 
    // used in the reflux operation
    flux_reg.resize(nlevs_max+1);
}

AmrCoreAdv::~AmrCoreAdv ()
{

}

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    BL_PROFILE("AmrCoreAdv::Evolve()");
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration);

	cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

	// sync up time
	for (int lev = 0; lev <= finest_level; ++lev) {
	    t_new[lev] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();
	}

	if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
	WritePlotFile();
    }
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    const Real time = 0.0;
    InitFromScratch(time);
    AverageDown();

    if (plot_int > 0) {
        WritePlotFile();
    }
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1]->nComp();
    const int nghost = phi_new[lev-1]->nGrow();
    
    MFInfo mfinfo;
#ifdef CUDA
    mfinfo.SetDevice(true);
#else
    mfinfo.SetDevice(false);
#endif
    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost, mfinfo));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost, mfinfo));


    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and 
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
		     const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev]->nComp();
    const int nghost = phi_new[lev]->nGrow();

    MFInfo mfinfo;
#ifdef CUDA
    mfinfo.SetDevice(true);
#else
    mfinfo.SetDevice(false);
#endif
#if __cplusplus >= 201402L
    auto new_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost, mfinfo);
    auto old_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost, mfinfo);
#else
    std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, ncomp, nghost, mfinfo));
    std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, ncomp, nghost, mfinfo));
#endif

    FillPatch(lev, time, *new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }    
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    phi_new[lev].reset(nullptr);
    phi_old[lev].reset(nullptr);
    flux_reg[lev].reset(nullptr);
}

// initialize data using a fortran routine to compute initial state
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				      const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 0;

    MFInfo mfinfo;
#ifdef CUDA
    mfinfo.SetDevice(true);
#else
    mfinfo.SetDevice(false);
#endif
    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost, mfinfo));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost, mfinfo));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

	initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
		 BL_TO_FORTRAN_3D(state[mfi]), ZFILL(dx),
		 ZFILL(prob_lo));
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Array<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
	first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
	ParmParse pp("adv");
	int n = pp.countval("phierr");
	if (n > 0) {
	    pp.getarr("phierr", phierr, 0, n);
	}
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = *phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
	
	for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box& tilebox  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    
	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
            // set itags initially to 'untagged' everywhere
            // we define itags over the tilebox region
	    tagfab.get_itags(itags, tilebox);
	    
            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebox.loVect();
	    const int*  thi     = tilebox.hiVect();

            // tag cells for refinement
	    state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
			BL_TO_FORTRAN_3D(state[mfi]),
			&tagval, &clearval, 
			ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()), 
			ZFILL(dx), ZFILL(prob_lo), &time, &phierr[lev]);
	    //
	    // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
	    //
	    tagfab.tags_and_untags(itags, tilebox);
	}
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("regrid_int", regrid_int);
	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);
    }

    {
	ParmParse pp("adv");
	
	pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(*phi_new[lev+1], *phi_new[lev],
			     geom[lev+1], geom[lev],
			     0, phi_new[lev]->nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(*phi_new[crse_lev+1], *phi_new[crse_lev],
			 geom[crse_lev+1], geom[crse_lev],
			 0, phi_new[crse_lev]->nComp(), refRatio(crse_lev));
}

// compute the number of cells at a level
long
AmrCoreAdv::CountCells (int lev)
{
    const int N = grids[lev].size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; ++i)
    {
        cnt += grids[lev][i].numPts();
    }

    return cnt;
}

// compute a new multifab by copying in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (my_verbose)
        std::cout << "AmrCoreAdv::FillPatch with lev = " << lev << std::endl;
    if (lev == 0)
    {
	Array<MultiFab*> smf;
	Array<Real> stime;
	GetData(0, time, smf, stime);

	AmrCoreAdvPhysBC physbc;
	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
				     geom[lev], physbc);
    }
    else
    {
	Array<MultiFab*> cmf, fmf;
	Array<Real> ctime, ftime;
	GetData(lev-1, time, cmf, ctime);
	GetData(lev  , time, fmf, ftime);

	AmrCoreAdvPhysBC cphysbc, fphysbc;
	Interpolater* mapper = &cell_cons_interp;

	int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
	int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
	Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

	amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
				   0, icomp, ncomp, geom[lev-1], geom[lev],
				   cphysbc, fphysbc, refRatio(lev-1),
				   mapper, bcs);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Array<MultiFab*> cmf;
    Array<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    
    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    AmrCoreAdvPhysBC cphysbc, fphysbc;
    Interpolater* mapper = &cell_cons_interp;
    
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, fphysbc, refRatio(lev-1),
				 mapper, bcs);
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrCoreAdv::GetData (int lev, Real time, Array<MultiFab*>& data, Array<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
	data.push_back(phi_new[lev].get());
	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
	data.push_back(phi_old[lev].get());
	datatime.push_back(t_old[lev]);
    }
    else
    {
	data.push_back(phi_old[lev].get());
	data.push_back(phi_new[lev].get());
	datatime.push_back(t_old[lev]);
	datatime.push_back(t_new[lev]);
    }
}


// advance a level by dt
// includes a recursive call for finer levels
void
AmrCoreAdv::timeStep (int lev, Real time, int iteration)
{
    BL_PROFILE("AmrCoreAdv::timeStep()");
    if (regrid_int > 0)  // We may need to regrid
    {

        BL_PROFILE("AmrCoreAdv::regrid()");
        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Array<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if 
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) 
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
		int old_finest = finest_level; 
                if (Verbose()) {
                    amrex::Print() << "Regrid based on level: " << lev << std::endl;
                }
		regrid(lev, time);

                // mark that we have regridded this level already
		for (int k = lev; k <= finest_level; ++k) {
		    last_regrid_step[k] = istep[k];
		}

                // if there are newly created levels, set the time step
		for (int k = old_finest+1; k <= finest_level; ++k) {
		    dt[k] = dt[k-1] / MaxRefRatio(k-1);
		}
	    }
	}
    }

    if (Verbose()) {
	amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
	amrex::Print() << "ADVANCE with dt = " << dt[lev] << std::endl;
    }

    // advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose())
    {
	amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
	for (int i = 1; i <= nsubsteps[lev+1]; ++i)
	{
	    timeStep(lev+1, time+(i-1)*dt[lev+1], i);
	}

	if (do_reflux)
	{
            // update lev based on coarse-fine flux mismatch
	    flux_reg[lev+1]->Reflux(*phi_new[lev], 1.0, 0, 0, phi_new[lev]->nComp(), geom[lev]);
	}

	AverageDownTo(lev); // average lev+1 down to lev
    }
    
}

// advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::Advance (int lev, Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("AmrCoreAdv::Advance()"); 
    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt;

    MultiFab& S_new = *phi_new[lev];

    if (my_verbose)
        std::cout << "Level: " << lev << " has " << S_new.size() << " FABs" << std::endl;
    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux)
    {
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    BoxArray ba = grids[lev];
	    ba.surroundingNodes(i);
            MFInfo mfinfo;
#ifdef CUDA
            mfinfo.SetDevice(true);
#else
            mfinfo.SetDevice(false);
#endif
	    fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0, mfinfo);
	}
    }

    MFInfo mfinfo;
#ifdef CUDA
    mfinfo.SetDevice(true);
#else
    mfinfo.SetDevice(false);
#endif
    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow, mfinfo);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef CUDA
    // std::vector<intptr_t> m_mem_tags;
    // m_mem_tags.reserve(n_fabs);
    // TODO: find a better solution to this
    intptr_t* m_mem_tags = new intptr_t[100000];
    std::array<std::vector<FArrayBox*>, BL_SPACEDIM> m_copy_dst;
    std::array<std::vector<Box>, BL_SPACEDIM> m_copy_box;
#endif

    std::array<std::vector<FArrayBox>, BL_SPACEDIM> flux_fabs;
    std::array<std::vector<FArrayBox>, BL_SPACEDIM> uface_fabs;

    BL_PROFILE_VAR("AmrCoreAdv::Advance()::advect_group_all", advect_group_all);
    BL_PROFILE_VAR("AmrCoreAdv::Advance()::advect_group_cpu", advect_group_cpu);
#ifdef _OPENMP
#pragma omp parallel 
#endif
    {
	FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

#ifdef CUDA
        // MFIter, tiling = true, use_device = true
        // Only when use_device = true, OMP thread 0 will get a MFIter
        // that has jobs for GPU
	// for (MFIter mfi(S_new, true, true); mfi.isValid(); ++mfi)
	for (MFIter mfi(S_new, true, true); mfi.isValid(); ++mfi)
#else
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
#endif
	{
	    const Box& bx = mfi.tilebox();

	    const FArrayBox& statein = Sborder[mfi];
	    FArrayBox& stateout      =   S_new[mfi];

            // TODO: write ifdef for these
            if ( omp_get_thread_num() == 0 ) { // only thread 0 talks to GPU
            // if ( false ) { // only thread 0 talks to GPU
                // int idx = mfi.LocalIndex();
                int idx = mfi.tileIndex();
                int dev_id = statein.deviceID();
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    const Box& bxtmp = amrex::surroundingNodes(bx,i);
                    flux[i].usePinnedMemory(true);
                    flux[i].needDeviceCopy(true);
                    flux[i].setDevice(0);
                    flux[i].resize(bxtmp,S_new.nComp());

                    uface[i].usePinnedMemory(true);
                    uface[i].needDeviceCopy(true);
                    uface[i].setDevice(0);
                    uface[i].resize(amrex::grow(bxtmp,1),1);
                }

                // compute velocities on faces (prescribed function of space and time)
                get_face_velocity(lev, ctr_time,
                                  AMREX_D_DECL(BL_TO_FORTRAN_DEVICE((uface[0])),
                                               BL_TO_FORTRAN_DEVICE((uface[1])),
                                               BL_TO_FORTRAN_DEVICE((uface[2]))),
                                  dx, prob_lo
                                  , idx, dev_id, idx);
                statein.toDevice(idx);
                advect(time, bx.loVect(), bx.hiVect(),
                       BL_TO_FORTRAN_3D_DEVICE(statein), 
                       BL_TO_FORTRAN_3D_DEVICE(stateout),
                       AMREX_D_DECL(BL_TO_FORTRAN_3D_DEVICE(uface[0]),
                                    BL_TO_FORTRAN_3D_DEVICE(uface[1]),
                                    BL_TO_FORTRAN_3D_DEVICE(uface[2])),
                       AMREX_D_DECL(BL_TO_FORTRAN_3D_DEVICE(flux[0]), 
                                    BL_TO_FORTRAN_3D_DEVICE(flux[1]), 
                                    BL_TO_FORTRAN_3D_DEVICE(flux[2])), 
                       dx, dt, idx, dev_id, idx);
                stateout.toHost(idx);
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    flux[i].toHost(idx);
                }
                if (do_reflux) {
                    for (int i = 0; i < BL_SPACEDIM ; i++) {
                        m_copy_dst[i].push_back(&(fluxes[i][mfi]));
                        m_copy_box[i].push_back(mfi.nodaltilebox(i));
                    }
                }
                // add callback function to CUDA stream associated with idx
                cudaStream_t pStream;
                get_stream(&idx, &pStream, &dev_id);
                m_mem_tags[idx] = idx;
                cudaStreamAddCallback(pStream, cudaCallback_release_gpu, (void*) &(m_mem_tags[idx]), 0);

                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    flux_fabs[i].push_back(std::move(flux[i]));
                    uface_fabs[i].push_back(std::move(uface[i]));
                }
                if (do_reflux) {
                    for (int i = 0; i < BL_SPACEDIM ; i++) {
                    }
                }
            } else {
            // BL_PROFILE_CUDA_GROUP("advect_group_cpu_exclusive", omp_get_thread_num());
            std::cout << "get_face_velocity_host" << std::endl;
                // Allocate fabs for fluxes and Godunov velocities.
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    const Box& bxtmp = amrex::surroundingNodes(bx,i);
                    flux[i].resize(bxtmp,S_new.nComp());
                    uface[i].resize(amrex::grow(bxtmp,1),1);
                }
                get_face_velocity_host(lev, ctr_time,
                                  AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                         BL_TO_FORTRAN(uface[1]),
                                         BL_TO_FORTRAN(uface[2])),
                                  dx, prob_lo);

                // compute new state (stateout) and fluxes.
                advect_host(time, bx.loVect(), bx.hiVect(),
                       BL_TO_FORTRAN_3D(statein), 
                       BL_TO_FORTRAN_3D(stateout),
                       AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
                              BL_TO_FORTRAN_3D(uface[1]),
                              BL_TO_FORTRAN_3D(uface[2])),
                       AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
                              BL_TO_FORTRAN_3D(flux[1]), 
                              BL_TO_FORTRAN_3D(flux[2])), 
                       dx, dt);
                if (do_reflux) {
                    for (int i = 0; i < BL_SPACEDIM ; i++) {
                        fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));	  
                    }
                }
            }

	}
    }
    BL_PROFILE_VAR_STOP(advect_group_cpu);
#ifdef CUDA
    // TODO: put this in destructor
    // also put destruction of m_mem_tags in destructor
    // In this way, we can create m_mem_tags after MFIter is constrctued,
    // by when we already know number of all tiles so we know how much 
    // memory we need to allocate for m_mem_tags
    int n_dev = ParallelDescriptor::get_num_devices_used();
    for (int i = 0; i < n_dev; ++i) {
        checkCudaErrors(cudaSetDevice(i));
        checkCudaErrors(cudaDeviceSynchronize());
    }
    if (do_reflux) {
        BL_PROFILE("AmrCoreAdv::Advance()::do_reflux");
        int N = flux_fabs[0].size();
        if ( N > 0 ) {
            for (int nf = 0; nf < N; ++nf) {
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    m_copy_dst[i][nf]->copy(flux_fabs[i][nf], m_copy_box[i][nf]);
                }
            }
        }
    }
#endif
    BL_PROFILE_VAR_STOP(advect_group_all);


    // increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // In this example we are solving phi_t = -div(+F)
    // The fluxes contain, e.g., F_{i+1/2,j} = (phi*u)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (do_reflux) { 
	if (flux_reg[lev+1]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
	        // update the lev+1/lev flux register (index lev+1)   
	        flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
	    }	    
	}
	if (flux_reg[lev]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
	        // update the lev/lev-1 flux register (index lev) 
		flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
	    }
	}
    }
}

// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
    Array<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
	dt_tmp[lev] = EstTimeStep(lev, true);
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
	dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
	n_factor *= nsubsteps[lev];
	dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
	dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
	dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, bool local) const
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];
    const MultiFab& S_new = *phi_new[lev];
    if (my_verbose) {
        std::cout << "In AmrCoreAdv::EstTimeStep()" << std::endl;
        std::cout << "Lev: " << lev << " has " << S_new.size() << " FABs." << std::endl;
    }

#ifdef CUDA
    // store tag for each fab here
    int n_fabs = S_new.size();
    intptr_t* m_mem_tags = new intptr_t[n_fabs];
#endif

    MultiFab ufaces[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        BoxArray ba = grids[lev];
        ba.surroundingNodes(i);
        MFInfo mfinfo;
#ifdef CUDA
        mfinfo.SetDevice(true);
#else
        mfinfo.SetDevice(false);
#endif
        ba.grow(1);
        ufaces[i].define(ba, dmap[lev], S_new.nComp(), 0, mfinfo);
    }


#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	FArrayBox uface[BL_SPACEDIM];


#ifdef CUDA
        // MFIter, tiling = false, use_device = true
	// for (MFIter mfi(S_new, false, true); mfi.isValid(); ++mfi)
	for (MFIter mfi(S_new, false, true); mfi.isValid(); ++mfi)
#else
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
#endif
	{
#ifdef CUDA
            int idx = mfi.LocalIndex();
            int dev_id = ufaces[0][mfi].deviceID();
            m_mem_tags[idx] = (intptr_t) &(S_new[mfi]);
#endif

            // TODO: write ifdef for these
            if ( omp_get_thread_num() == 0 ) { // only thread 0 talks to GPU
                get_face_velocity(lev, cur_time,
                                  AMREX_D_DECL(BL_TO_FORTRAN_DEVICE(ufaces[0][mfi]),
                                         BL_TO_FORTRAN_DEVICE(ufaces[1][mfi]),
                                         BL_TO_FORTRAN_DEVICE(ufaces[2][mfi])),
                                  dx, prob_lo
                                  , idx, dev_id, m_mem_tags[idx]
                                  );
                // add callback function to CUDA stream associated with idx
                cudaStream_t pStream;
                get_stream(&idx, &pStream, &dev_id);
                cudaStreamAddCallback(pStream, cudaCallback_release_gpu, (void*) &m_mem_tags[idx], 0);
                dt_est = std::numeric_limits<Real>::max(); // will process this later
                if (my_verbose) {
#pragma omp critical
                    std::cout << "Thread: " << omp_get_thread_num() << " works on FAB: " << idx << ". (on GPU)" << std::endl;
                }
             } else {
//                 if (my_verbose) {
// #pragma omp critical
//                     std::cout << "Thread: " << omp_get_thread_num() << " works on FAB: " << idx << std::endl;
//                 }
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    const Box& bx = mfi.nodaltilebox(i);
                    uface[i].resize(bx,1);
                }
                get_face_velocity_host(lev, cur_time,
                                  AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                         BL_TO_FORTRAN(uface[1]),
                                         BL_TO_FORTRAN(uface[2])),
                                  dx, prob_lo);
                for (int i = 0; i < BL_SPACEDIM; ++i) {
                    Real umax = uface[i].norm(0);
                    if (umax > 1.e-100) {
                        dt_est = std::min(dt_est, dx[i] / umax);
                    }
                }
             }

	}
    }
#ifdef CUDA
    // synchronize all devices
    int n_dev = ParallelDescriptor::get_num_devices_used();
    for (int i = 0; i < n_dev; ++i) {
        checkCudaErrors(cudaSetDevice(i));
        checkCudaErrors(cudaDeviceSynchronize());
    }
    Real* umax = (Real*) malloc(1 * sizeof(Real));

// TODO: right now if you don't use Openmp here, you can't filter out 
// FABs that needs these processing
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, false, true); mfi.isValid(); ++mfi) {
        // TODO: write ifdef for these
        // compute norm of this fab on GPU
        if ( omp_get_thread_num() == 0 ) { // Only FABs processed by GPU need the following
            if (my_verbose) {
#pragma omp critical
                std::cout << "Thread: " << omp_get_thread_num() << " call cublasIdamax on FAB: " << mfi.LocalIndex() << ". (on GPU)" << std::endl;
            }
            for (int i = 0; i < BL_SPACEDIM; ++i) {
                int max_id = -1;
                // assume we only have one component in the FAB now
                // assume we use only one GPU
                BL_ASSERT(ParallelDescriptor::get_num_devices_used() == 1);
                cublasIdamax(cublasHandles[0], ufaces[i][mfi].nPts(), ufaces[i][mfi].devicePtr(), 1, &max_id);
                checkCudaErrors(cudaMemcpy(umax, &((ufaces[i][mfi].devicePtr())[max_id]), 1*sizeof(Real),cudaMemcpyDeviceToHost));
                if (*umax > 1.e-100) {
                    dt_est = std::min(dt_est, dx[i] / (*umax));
                }
            }
        }
    }
    delete[] m_mem_tags;
#endif

    if (!local) {
	ParallelDescriptor::ReduceRealMin(dt_est);
    }

    dt_est *= cfl;

    return dt_est;
}

// get plotfile name
std::string
AmrCoreAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Array<const MultiFab*>
AmrCoreAdv::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(phi_new[i].get());
    }
    return r;
}

// set plotfile variable names
Array<std::string>
AmrCoreAdv::PlotFileVarNames () const
{
    return {"phi"};
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				    Geom(), t_new[0], istep, refRatio());
}
