#!/bin/bash

#=========== Global directives ===========
# @ environment = COPY_ALL
# @ job_name = compile_wrf
# @ output = $(job_name).$(jobid)
# @ error  = $(job_name).$(jobid)
# @ job_type = serial
# @ parallel_threads = 8
# @ requirements = ( Feature == "prepost" )
# @ as_limit = 20.0gb
# @ wall_clock_limit = 03:00:00
# @ queue

set -x

./make_WRF_compil
