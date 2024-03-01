
#########################   from common_definitions.sh  ########################
#-------------------------------------------------------------------------------
#  Default values
#-------------------------------------------------------------------------------
 export MACHINE_STOCKAGE=""
#export MACHINE_STOCKAGE="ssh gaya"

#-------------------------------------------------------------------------------
#  Common functions
#-------------------------------------------------------------------------------
function lnfile  # (if the file is unreachable we stop: EXIT)
{
  [ "$#" -ne 2 ] && printf "\n\nlnfile $* \n\n wrong (needs 2 files), we stop...\n\n" && exit 1
  [ ! -f $1 ] && printf "\n\nlnfile: file $1 does not exist, we stop...\n\n" && exit 1
  ln -s $1 $2
  echo "ln -s $1 -> $2 ok"
}

function cpfile  # (if the file is unreachable we stop: EXIT)
{
  [ "$#" -ne 2 ] && printf "\n\ncpfile $* \n\n wrong (needs 2 files), we stop...\n\n" && exit 1
  [ ! -f $1 ] && printf "\n\ncpfile: file $1 does not exist, we stop...\n\n" && exit 1
  cp -f $1 $2
  echo "cp $1 -> $2 ok"
}

function cpfile2  # (even if the file is unreachable we continue)
{
  [ "$#" -ne 2 ] && printf "\n\ncpfile2 $* \n\n wrong (needs 2 files), we stop...\n\n" && exit 1
  if [ -f $1 ]; then
    cp -f $1 $2
    echo "cp $1 -> $2 ok"
  else
    printf "cpfile2: file $1 does not exist, we continue...\n"
  fi
}

function mvfile   # (if the file is unreachable we stop: EXIT)
{
  [ "$#" -ne 2 ] && printf "\n\nmvfile $* \n\n wrong (needs 2 files), we stop...\n\n" && exit 1
  [ ! -f $1 ] && printf "\n\nmvfile: file $1 is not existing, we stop...\n\n" && exit 1
  mv -f $1 $2
  echo "mv $1 -> $2 ok"
}

function mvfile2   # (even if the file is unreachable we continue)
{
  [ "$#" -ne 2 ] && printf "\n\nmvfile2 $* \n\n wrong (needs 2 files), we stop...\n\n" && exit 1
  if [ -f $1 ]; then
    mv -f $1 $2
    echo "mv $1 -> $2 ok"
  else
    printf "mvfile2: file $1 does not exist, we continue...\n"
  fi
}

