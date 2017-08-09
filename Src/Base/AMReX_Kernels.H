#ifndef BL_KERNELS_H_ 
#define BL_KERNELS_H_
#include <AMReX_FabArray.H>
#define ARRAY_2D(PHI, LO_X, LO_Y, HI_X, HI_Y, I, J) PHI[(J-LO_Y)*(HI_X-LO_X+1)+I-LO_X]
    // template <class FAB> class FabArray;
    // template <class FAB> using TagType = typename FabArray<FAB>::CopyComTagDevice;
    // template<class FAB>
    // void fill_bc_d2d(TagType<FAB>* tags, int n_tags, int scomp, int ncomp) {return;} 


    // template<class FAB> class FabArray;
    // // __global__
    // template<class FAB>
    // void fill_bc_d2d(typename FabArray<FAB>::CopyComTagDevice* tags, int n_tags, int scomp, int ncomp) {
    //     return;
    // }
        // const int i_tag = blockIdx.y;
        // const int tid = threadIdx.x + blockDim.x * blockIdx.x;
        // const typename FabArray<FAB>::CopyComTagDevice& tag = tags[i_tag];
        // // unpacking
        // amrex::Real* src = tag.sfab_data_pt_d;
        // amrex::Real* dst = tag.dfab_data_pt_d;
        // const int sfab_lox = tag.sfab_box.loVect()[0];
        // const int sfab_loy = tag.sfab_box.loVect()[1];
        // const int sfab_hix = tag.sfab_box.hiVect()[0];
        // const int sfab_hiy = tag.sfab_box.hiVect()[1];

        // const int dfab_lox = tag.dfab_box.loVect()[0];
        // const int dfab_loy = tag.dfab_box.loVect()[1];
        // const int dfab_hix = tag.dfab_box.hiVect()[0];
        // const int dfab_hiy = tag.dfab_box.hiVect()[1];

        // const int sbox_lox = tag.sbox.loVect()[0];
        // const int sbox_loy = tag.sbox.loVect()[1];
        // const int sbox_hix = tag.sbox.hiVect()[0];
        // const int sbox_hiy = tag.sbox.hiVect()[1];

        // const int dbox_lox = tag.dbox.loVect()[0];
        // const int dbox_loy = tag.dbox.loVect()[2];
        // const int dbox_hix = tag.dbox.hiVect()[0];
        // const int dbox_hiy = tag.dbox.hiVect()[1];

        // //TODO: assume sbox == dbox for now
        // const int boxSize_x = sbox_hix - sbox_lox + 1;
        // const int x_local = tid/boxSize_x + 1;
        // const int y_local = tid%boxSize_x;
        // const int x_global = x_local + sbox_lox;
        // const int y_global = y_local + sbox_loy;
        // if (x_global <= sbox_hix || y_global <= sbox_hiy ) {
        //     ARRAY_2D(dst,dfab_lox,dfab_loy,dfab_hix,dfab_hiy,x_global,y_global) = 
        //         ARRAY_2D(src,sfab_lox,sfab_loy,sfab_hix,sfab_hiy,x_global,y_global);  
        // }
    // }

#endif
