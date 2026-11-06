#!/usr/bin/env bash
# Report an LSF job's state to Snakemake: running, success or failed.
jobid=$1
stat=$(bjobs -noheader -o stat "$jobid" 2>/dev/null | tail -n1)
case "$stat" in
    DONE) echo success ;;
    EXIT) echo failed ;;
    PEND|RUN|PSUSP|USUSP|SSUSP|WAIT|PROV) echo running ;;
    *)  # no longer listed by bjobs: look in the job history
        if bhist -l "$jobid" 2>/dev/null | grep -q "Done successfully"; then echo success
        elif bhist -l "$jobid" 2>/dev/null | grep -qE "Exited|Completed <exit>"; then echo failed
        else echo running
        fi ;;
esac
