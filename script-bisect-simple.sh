#!/bin/bash -fx

# faf26e526c7f7e38f819b2d93ad161eb01f1413f



current_last_commit=`git log -1 --format="%H"`
echo $current_last_commit
git log -1


current_commit_date=`git log -1 --pretty=format:"%cd" --date=short`
echo $current_commit_date




#commit_make_mechanisms='87dd17ad462c4b36f6544d92c14ea79147d49e2c'
#echo $commit_make_mechanisms
# if (git rev-list $current_last_commit| grep $(git rev-parse $commit_make_mechanisms))
# then
#     echo "We are after make mechanisms as a component"
#     cd /scratch/Vincent/Build/siconos
#     \rm -rf bisect
#     mkdir bisect
#     cd bisect
#     ln -sf ~/build_options_new.cmake ~/build_options.cmake
#     make -f  ~/build.Makefile siconos_debug
# else
#     echo "We are before make mechanisms as a component"
#     cd /scratch/Vincent/Build/siconos
#     \rm -rf bisect
#     mkdir bisect
#     cd bisect
#     ln -sf ~/build_options_old.cmake ~/build_options.cmake
#     make -f  ~/build.Makefile siconos_debug
# fi

# make -j8 install || exit 125
cd /Users/acary/siconos.old/kernel/src
#sed -e "s/\#include\ <boost\/timer.hpp>//g" SiconosKernel.hpp > SiconosKernel.hpp.wait
#sed -e "s/\#include\ <boost\/progress.hpp>//g" SiconosKernel.hpp.wait > SiconosKernel.hpp


cd /scratch/Vincent/Build/siconos
\rm -rf bisect
mkdir bisect
cd bisect
#make -f  ~/build.Makefile siconos_debug; make -j8 install || exit 125
make -f ~/build.Makefile  siconos_clean;
make -f ~/build.Makefile  siconos_old_kernel_debug; make -j8 install 2> /dev/null || exit 125

#cd /Users/acary/siconos.old/kernel/src
#git ch SiconosKernel.hpp

cd /scratch/vincent/siconos-tutorials.old/examples/mechanics/JointsTests
commit_tutorials=`git rev-list -1 --before=$current_commit_date master --format=medium`
echo $commit_tutorials
git checkout `git rev-list -1 --before=$current_commit_date master` -b `git rev-list -1 --before=$current_commit_date master`
git lg -1

\rm -rf .siconos NE_1DS_1Knee_MLCP
#siconos --clean
#siconos -c

siconos NE_1DS_1Knee_MLCP.cpp
status=$?

if [ $status -eq 0 ]
then
    echo $status
    exit $status
else
    git log -1
    echo $status
    exit 1
fi
