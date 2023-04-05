
#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include "AmrMHD.H"

using namespace amrex;

amrex::LevelBld* getLevelBld ();

int
main (int   argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    Real dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    Real stop_time;
    bool test_mode;

    {
        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step", max_step);
        pp.query("strt_time", strt_time);
        pp.query("stop_time", stop_time);
        pp.query("test_mode", test_mode);

    }

    if (test_mode) {
        // State uL = {1.5, 0.5, 3.0, 4.0, 3.0, 2.2, 0.8};
        // State uR = {4.0, 0.3, 7.0, 1.9, 2.5, 5.6, 0.5};
        State fluxVec;
        // HLLCFlux(uL, uR, fluxVec, 0);
        // for (int i = 0; i < stateVars; i++) {
        //     amrex::Print() << fluxVec[i] << " ";
        // }
        // amrex::Print() << std::endl;
        State wL = {1., 0., 0., 0., 1., 0.75, 1., 0.};
        State wR = {0.125, 0., 0., 0., 0.1, 0.75, -1, 0.};
        State uL, uR;
        conservative(wL, uL);
        conservative(wR, uR);

        // HLLCFlux(uL, uL, fluxVec, 0);
        amrex::Print() << "uL = ";
        for (int i = 0; i < stateVars; i++) {
            amrex::Print() << uL[i] << " ";
        }
        amrex::Print() << std::endl;
        amrex::Print() << "uR = ";
        for (int i = 0; i < stateVars; i++) {
            amrex::Print() << uR[i] << " ";
        }
        amrex::Print() << std::endl;
    }
    else {
        if (strt_time < 0.0) {
            amrex::Abort("MUST SPECIFY a non-negative strt_time");
        }

        if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
        }

        {
            Amr amr(getLevelBld());

        amr.init(strt_time,stop_time);

        while ( amr.okToContinue() &&
            (amr.levelSteps(0) < max_step || max_step < 0) &&
            (amr.cumTime() < stop_time || stop_time < 0.0) )

        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
            //
            amr.coarseTimeStep(stop_time);
        }

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }

        }

        Real dRunTime2 = amrex::second() - dRunTime1;

        ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "Run time = " << dRunTime2 << std::endl;
    }

    amrex::Finalize();

    return 0;
}
