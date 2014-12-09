#!/bin/bash

# Build release version of Mantid
cd $HOME/repositories/mantidproject/mantid/
git stash
git checkout master
git fetch -p
git pull --rebase
cmakeMantid.py --builddir release --ncore 1 --popmsg no
cmakeMantid.py --builddir release --ncore 1 --target doc --popmsg no
git checkout -
git stash apply
