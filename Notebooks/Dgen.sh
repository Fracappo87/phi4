#!/bin/bash

#INSTRUCTION: function that gives instructions about input parameters.


INSTRUCTION(){
  cat <<EOF  
INSTRUCTIONS: $0 [-h] [-D <max number derivative>]  

OPTIONS:
   -h      Show this message
   -D      Max number of derivative
The order all this options are given is not important.

EOF
}





#Initial definition of parameters that will be used during the script


D=
codefile="derivative.m"
codefilejac="derivative_jac.m"
HOMEDIR=$HOME

THISSCRIPTb=`basename $0`

#GETOPT: start to accept input flags with their related values: depending on which flags are stord in OPTIONS, the script will do different things



while getopts "hD:" OPTION
do
    case $OPTION in
        h)
            INSTRUCTION
            exit 1
            ;;
        D)
            D=$OPTARG
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



if [[ -z $D ]] 
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
ext="macro.h"

./$codefile $D Phiboundary Phi > $ext 

mv $ext $INCLUDE/

ext="macro_jacobian.h"

./$codefilejac $D Jacobian > $ext 

mv $ext $INCLUDE/

echo ""
echo -e "$THISSCRIPT: evaluation terminated \n"
echo ""

    

exit $?

