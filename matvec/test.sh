#!/bin/bash
sci_notation_regex='^[0-9]+([.][0-9]+)?(e[0-9]+|e-[0-9]+)?$'

function test_time {
    # compare 
    if [[ ! $1 =~ $sci_notation_regex ]] ; 
    then
        echo ERROR: time is not on stderr or not formatted properly
        echo
        rm .time
        exit 1
    fi
    # delete tmp file 
    rm .time
}

SUCCESS_FILE=.passed_mpi_matmul
# remove success file
if [ -e ${SUCCESS_FILE} ] ; 
then
    rm ${SUCCESS_FILE} 
fi


. ../params.sh

NS="720 2880"
ITERS="1 3 4"
NPP="4 9"
for NP in ${NPP} ;
do
   for N in ${NS} ;
   do
      for ITER in ${ITERS} ;
      do

         ANSW=$(mpirun --mca btl_base_warn_component_unused 0 ${MPIRUN_PARAMS} -np ${NP} ./mpi_matmul ${N} ${ITER} 2> .time < /dev/null)

	 process_time_file .time
	 
	 #The test use to read: I believe this is wrong because it reruns the code sequentially (Erik:2019/04/01)
#         if [ -z "$(./mpi_matmul ${N} ${ITER} 2> .time)" ] ;

	 if [ -z "${ANSW}" ] ;
         then
            test_time $(cat .time) 
         else
            echo FAILED: "mpirun --mca btl_base_warn_component_unused 0 ${MPIRUN_PARAMS} -np ${NP} ./mpi_matmul ${N} ${ITER}"
            exit 1
         fi
      
       # "progess bar"
       echo -n "|"
      
      done
   done
done


touch ${SUCCESS_FILE}
echo
