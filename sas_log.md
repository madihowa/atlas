# Changes Made by SAS (06/18/21)

- moved all *.e* and *.o* into quanah_logs directory to clean up base directory
- new *.e* and *.o* will be generated in quanah_logs -> changed header of matador_job.sh
- added print out of run number and date in the Master.py so that we can see that info in the *.o* files and match it back to
the jobs that we ran on quanah
- added new aliases - clc, qj, ql - to go to ClusterCalbration, quanah_jobs and quanah_log directory straight
- added new aliases -sq and sa- to replace "squeue -u madihowa" and "sacct -j" respectively
- trying debugging idea 1
    - worked with logged on
    - didnt work when not logged on
- changed .ssh/config on both client and server didnt work
-[] need to email quanah about this 
- trying setting up XDG_RUNTIME_PATH like quanah people said
    - updated run.sh for this
    - fixed issue of XDG_RUN_PATH 
    - X error still persistent
    - trying -X only
        - X error still persistent
    - trying -Y only
        - X error still persistent
- trying with no .Xauthority
    - X error still persistent
- trying with QT_QPA_PLATFORM='offscreen'
    - problem solved!
- running a 100 epoch run! 
    - folder name: Results_2021-06-18_run_11/
- root now loads on startup -> added to .bashrc
- added last line to PlotHisto.C to terminate root after generating pics
- QUANAH IS NOW READY FOR BATCH JOBS FOR THIS SOFTWARE



# Debugging Ideas:
-[x] test running a job while logged onto quanah vs not logged on for epoch = 5
-[x] noticed different in error messages of display variable when quanah accessed through MIH laptop and SAS laptop
    - may need to set $DISPLAY to the location of error message
-[x] try setting XDG_RUNTIME_PATH to some actual path like quanah people asked with right permissions
-[x] try -X only and -Y only
-[x] delete .Xauthority
-[x] set QT_QPA_PLATFORM='offscreen'
