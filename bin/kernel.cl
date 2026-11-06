//
// Created by lishengting on 2019/9/29.
//
//#pragma OPENCL EXTENSION cl_nv_printf : enable

__kernel void sw (__global GI *gin, __global GO *gout) {
    int ii = get_global_id(0);
    //printf("<=%d=>\n",x);
    globalAlignAffine_seq(
            gin[ii].xseq,
            gin[ii].yseq,
            &gin[ii].aqmp,
            gin[ii].min_mlen,
            gin[ii].max_mlen,
            gout[ii].align_rets
            );
}

#if defined(TWO_TRACE) && defined(CL_TRACE)
__kernel void trace (__global TI1 *tin1, __global TI2 *tin2, __global TO *tout) {
    int ii = get_global_id(0);
    //printf("<=%d=>\n",x);
    int jj = tin2[ii].idx;
    traceback(
            tin1[jj].xseq,
            tin1[jj].yseq,
            &tin2[ii].align_ret,
            &tout[ii].align_ret
            );
}
#endif