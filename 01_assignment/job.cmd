#!/bin/bash
#@ job_name = pos_ag
#@ job_type = parallel
#@ output = job.out
#@ error = job.err
#@ class = test
#@ node = 1
#@ total_tasks = 16
#@ energy_policy_tag = NONE
#@ island_count = 1
#@ wall_clock_limit = 01:00:00
#@ notification = always
#@ notify_user = grant.bartel@tum.de
#@ queue
./gccg bin data/pent.bin pent
echo "The job has finished!"