
#include "Nyx.H"
#include "Nyx_error_F.H"

using namespace amrex;

void
Nyx::error_setup()
{
    // The lines below define routines to be called to tag cells for error
    // estimation -- the arguments of each "add" call are:
    //   1. Name of variable (state variable or derived quantity) which will be
    //      passed into the Fortran subroutine.
    //   2. Number of ghost cells each array needs in each call to the Fortran
    //      subroutine
    //   3. Type of Fortran subroutine -- this determines the argument list of
    //      the Fortran subroutine. These types are pre-defined and are
    //      currently restricted to `ErrorRec::Standard` and
    //      `ErrorRec::UseAverage`.
    //   4. Name of Fortran subroutine.

    err_list.add("total_particle_count", 1, ErrorRec::Standard,
                 BL_FORT_PROC_CALL(TAG_PART_CNT_ERR, tag_part_cnt_err));
}

void
Nyx::manual_tags_placement (TagBoxArray&    tags,
                            const Vector<IntVect>& bf_lev)
{
}
