#!/bin/bash

# author: Duane Johnson
# email: duane.johnson@gmail.com
# date: 2008 Jun 12
# license: MIT
# 
# Based on discussion at http://kerneltrap.org/mailarchive/git/2007/11/12/406496

pushd . >/dev/null

cd $1

# Find base of git directory
while [ ! -d .git ] && [ ! `pwd` = "/" ]; do cd ..; done

# Show various information about this git directory
if [ -d .git ]; then
  echo "== Remote URL: `git remote -v`"

#  echo "== Remote Branches: "
#  git branch -r
#  echo

  echo "== current revision of HEAD: "
  git rev-list HEAD | wc -l
  echo

  echo "== Local Branches:"
  git branch
  echo

#  echo "== Configuration (.git/config)"
#  cat .git/config
#  echo

  echo "== Most Recent Commit"
  git log --max-count=1
  echo

  echo "Type 'git log' for more commits, or 'git show' for full commit details."
else
  echo "Not a git repository."
  if [ ${USE_ATM} == 1 ]; then
      if [ $1 == ${ATM} ]; then
         [ -f README ]  && printf "\n $(head -n 1 README) \n" || printf "\n Could not find version \n"
      fi 
  fi
  if [ ${USE_OCE} == 1 ]; then
      if [ $1 == ${OCE} ]; then
         [ -f ../readme_version_croco.txt ]  && printf "\n $(head -n 3 readme_version_croco.txt) \n" || printf "\n Could not find version \n"
      fi
  fi
  if [ ${USE_WAV} == 1 ]; then
      if [ $1 == ${WAV} ]; then
         [ -f README ]  && printf "\n $(head -n 1 README) \n" || printf "\n Could not find version \n"
      fi
  fi
fi

popd >/dev/null
