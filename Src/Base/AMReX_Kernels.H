#ifndef BL_KERNELS_H_ 
#define BL_KERNELS_H_
#include <AMReX_FabArrayBase.H>
#define ARRAY_2D(PHI, LO_X, LO_Y, HI_X, HI_Y, I, J) PHI[(J-LO_Y)*(HI_X-LO_X+1)+I-LO_X]
namespace amrex {
    // template <class FAB> class FabArray;
    // template <class FAB> using TagType = typename FabArray<FAB>::CopyComTagDevice;
    // template<class FAB>
    // void fill_bc_d2d(TagType<FAB>* tags, int n_tags, int scomp, int ncomp) {return;} 
    //
    // TODO: always assume 2D and ncomp is 1 for now
    // class FabArrayBase;
    template<typename T>
    __global__
    void fill_bc_d2d(FabArrayBase::CopyComTagDevice* tags, int n_tags, int scomp, int ncomp) {
        const int i_tag = blockIdx.y;
        const int tid = threadIdx.x + blockDim.x * blockIdx.x;
        const FabArrayBase::CopyComTagDevice& tag = tags[i_tag];
        // unpacking
        // TODO: this should be more generic than Real
        T* src = (T*) tag.sfab_data_pt_d;
        T* dst = (T*) tag.dfab_data_pt_d;
        const int sfab_lox = tag.sfab_box.loVect()[0];
        const int sfab_loy = tag.sfab_box.loVect()[1];
        const int sfab_hix = tag.sfab_box.hiVect()[0];
        const int sfab_hiy = tag.sfab_box.hiVect()[1];

        const int dfab_lox = tag.dfab_box.loVect()[0];
        const int dfab_loy = tag.dfab_box.loVect()[1];
        const int dfab_hix = tag.dfab_box.hiVect()[0];
        const int dfab_hiy = tag.dfab_box.hiVect()[1];

        const int sbox_lox = tag.sbox.loVect()[0];
        const int sbox_loy = tag.sbox.loVect()[1];
        const int sbox_hix = tag.sbox.hiVect()[0];
        const int sbox_hiy = tag.sbox.hiVect()[1];

        const int dbox_lox = tag.dbox.loVect()[0];
        const int dbox_loy = tag.dbox.loVect()[1];
        const int dbox_hix = tag.dbox.hiVect()[0];
        const int dbox_hiy = tag.dbox.hiVect()[1];

        const int boxSize_x = dbox_hix - dbox_lox + 1;
        const int x_local = tid%boxSize_x;
        const int y_local = tid/boxSize_x;
        const int x_global = x_local + dbox_lox;
        const int y_global = y_local + dbox_loy;

        const int offx = sbox_lox - dbox_lox;
        const int offy = sbox_loy - dbox_loy;

        if (x_global <= dbox_hix && y_global <= dbox_hiy ) {
            ARRAY_2D(dst,dfab_lox,dfab_loy,dfab_hix,dfab_hiy,x_global,y_global) = 
                ARRAY_2D(src,sfab_lox,sfab_loy,sfab_hix,sfab_hiy,x_global+offx,y_global+offy);  
        }
    }

}

#endif
