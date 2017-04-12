#!/bin/bash

#INSTRUCTION: function that gives instructions about input parameters.


INSTRUCTION(){
  cat <<EOF  
INSTRUCTIONS: $0 [-h] [-D <max probe dimension>]  [-f <include file>]  [-p <probefile>]  [-v <varfile>] 

OPTIONS:
   -h      Show this message
   -D      Max dimension
   -f      Header file to be included in TWI_obs_engine.h
   -p      Probe file
   -v      Probe variation file
The order all this options are given is not important.

EOF
}





#Initial definition of parameters that will be used during the script


D=
f=
p=
v=
codefile="new_probe_gen.m"

HOMEDIR=$HOME

THISSCRIPTb=`basename $0`

#GETOPT: start to accept input flags with their related values: depending on which flags are stord in OPTIONS, the script will do different things



while getopts "hD:f:p:v:" OPTION
do
    case $OPTION in
        h)
            INSTRUCTION
            exit 1
            ;;
        D)
            D=$OPTARG
            ;;
        f)
            f=$OPTARG
            ;;
        p) 
            p=$OPTARG
            ;;
	v) 
            v=$OPTARG
            ;;
        ?)
            INSTRUCTION
            exit 1
            ;;
        esac
done




#======================================================================================================================================================================================================

#CHECKING: in this section the script checks that all the given options are correct

# ===> Input parameters check <===
# ===> Evaluation code existence <===
# ===> Creating files containing the final data <===



if [[ -z $D ]] || [[ -z $f ]] || [[ -z $p ]] || [[ -z $v ]]
then
     INSTRUCTION
     exit 1
fi


if [[ ! -f $codefile ]]
then
     echo -e ""
     echo -e "$THISSCRIPT: Missing the code that generates probes and their variations \n"
     echo -e ""
     exit 1
fi

INCLUDE=../Include
OBS=../Modules/Observables

./$codefile $D $f $p $v

mv ${f}.h $INCLUDE/
mv ${p}.h $INCLUDE/
mv ${v}.h $INCLUDE/
mv ${p}.c $OBS/
mv ${v}.c $OBS/
mv probevec $OBS/


echo ""
echo -e "$THISSCRIPT: evaluation terminated \n"
echo ""

    

exit $?

