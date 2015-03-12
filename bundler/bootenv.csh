#!/bin/tcsh

#Set up the bunder home directory
unset tmp && set tmp = `echo $_ | awk '{print $2 }' /dev/stdin` #this is really only so complicated cause this files needs to be sources, cant use $0
set tmp = `dirname $tmp`
unset BUNDLER_HOME && setenv BUNDLER_HOME `cd $tmp && pwd`
echo BUNDLER_HOME is $BUNDLER_HOME

#Set Default Number of cores for SIFT
setenv SIFT_LIMIT 4

#Set up loading the shared libraries, this is done by the runbundler scripts but we may want to run bundler without using those scripts
if (!($?LD_LIBRARY_PATH)) then
	setenv LD_LIBRARY_PATH "$BUNDLER_HOME/bin"
else
	setenv LD_LIBRARY_PATH "$LD_LIBRARY_PATH\:$BUNDLER_HOME/bin"
endif 

#Alias to run bundler from bundler git clone
alias RunBundler '$BUNDLER_HOME/RunBundler.sh'
alias RunBundler_profile '$BUNDLER_HOME/RunBundler_profile.sh'

alias top 'cd $BUNDLER_HOME' 
alias ds 'cd $BUNDLER_HOME/../data_sets'
alias res 'cd $BUNDLER_HOME/../results'

