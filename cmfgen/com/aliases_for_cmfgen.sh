#
# Put these aliases in you .cshrc file. If in your directory,
# you will need to replace hillier by you username
#
##############################################################
##############################################################
#
# cmfdist will be defined for general use (e.g. cd cmfdist)
# Change to you local distribution directory

# Altered 22-Jun-2022 - Adjusted to require less changing with updates
#                         and the use of GITHUB.

if ($machine =~ *'osiris'*)then
  set atomic  = ~hillier/cmfgen_atomic_data/atomic_9may23
  set ATOMIC   = ~hillier/cmfgen_atomic_data/atomic_9may23
  set sn_atomic  = ~hillier/sn_atomic 
  set SN_ATOMIC  = ~hillier/sn_atomic 
  if ($machine =~ *'wnode04'*)then
    setenv cmfdist  ~hillier/wnode04/cur_cmf
    setenv CMFDIST  ~hillier/wnode04/cur_cmf
  else
    setenv cmfdist  ~hillier/cur_cmf
    setenv CMFDIST  ~hillier/cur_cmf
  endif
else if ($machine =~ "ibis.phyast.pitt.edu" )then
  setenv cmfdist  ~jdh/Dropbox/Ibis/cur_cmf
  setenv CMFDIST  ~jdh/Dropbox/Ibis/cur_cmf
  set atomic  = ~jdh/atomic 
  set ATOMIC  = ~jdh/atomic 
  set sn_atomic  = ~jdh/sn_atomic 
  set SN_ATOMIC  = ~jdh/sn_atomic 
else if ($machine =~ *'Emu'* )then
  setenv cmfdist  ~jdh/cur_cmf
  setenv CMFDIST  ~jdh/cur_cmf
  set atomic  = ~jdh/atomic 
  set ATOMIC  = ~jdh/atomic 
  set sn_atomic  = ~jdh/sn_atomic 
  set SN_ATOMIC  = ~jdh/sn_atomic 
else
  setenv cmfdist  $home/cur_cmf
  setenv CMFDIST  $home/cur_cmf
  set atomic  = $home/atomic 
  set ATOMIC  = $home/atomic 
  set sn_atomic  = $home/sn_atomic 
  set SN_ATOMIC  = $home/sn_atomic 
endif

  set atomic  = ~hillier/cmfgen_atomic_data/atomic_9may23
  set ATOMIC   = ~hillier/cmfgen_atomic_data/atomic_9may23
  set sn_atomic  = ~hillier/cmfgen_atomic_data/sn_atomic 
  set SN_ATOMIC  = ~hillier/cmfgen_atomic_data/sn_atomic 
  setenv cmfdist  ~hillier/cur_cmf
  setenv CMFDIST  ~hillier/cur_cmf
  
#  echo 'Name of machine not recognized'
#  echo 'Update cmfdist/com/aliases_for_cmfgen.sh'


# Used to assign help fles for DISPGEN, PLT_SPEC etc
# Creates soft links.
# Must be executed inplotting directiry.

alias astxt $cmfdist/com/assign_txt_files.sh

#
# Executable shorthands

alias cmfgen         $cmfdist/exe/cmfgen.exe
alias cmf_flux       $cmfdist/exe/cmf_flux.exe
alias dispgen        $cmfdist/exe/dispgen.exe
alias plt_spec       $cmfdist/exe/plt_spec.exe
alias plt_scr        $cmfdist/exe/plt_scr.exe
alias plt_jh         $cmfdist/exe/plt_jh.exe

alias append_dc      $cmfdist/exe/append_dc.exe
alias rewrite_dc     $cmfdist/exe/rewrite_dc.exe
alias n_col_merge    $cmfdist/exe/n_col_merge.exe
alias n_pair_merge   $cmfdist/exe/n_pair_merge.exe

#
# Performs an NG acceleration operating on the SCRTEMP file

alias do_ng        $cmfdist/exe/do_ng.exe


#
# Used to delete uneccesary files (e.g. Scratch files) after Model has 
# completed. 

alias clean        $cmfdist/com/clean.sh

#
# Used to rename *OUT files to *_IN

alias out2in       $cmfdist/com/out2in.sh

#
# Used to name *.for files to *.f

alias for2f        $cmfdist/com/for_to_f.sh

#
# Used to name xxx.inc file to XXX.INC

alias inc2inc      $cmfdist/com/inc2inc.sh

##################################################
#
# Used to compare *.f files in one directory to those in another

alias for_dif        $cmfdist/com/for_dif.sh

#
# Used to compare *.INC files in one directory to those in another

alias INC_dif        $cmfdist/com/INC_dif.sh

#
# Used to copy files required for the generation of a new model from a 
# a directory containg a completed model.

alias cpmod       $cmfdist/com/cpmod.sh

# General aliases to clean a model directory. Be careful --- these might
# delete needed files if used incorrectly.

alias rmrrr     'rm -vf *PRRR'
alias rmin      'rm -vf *_IN'
alias dfort     'rm -vf fort.*'
alias dsve      'rm -vf *.sve'
alias dlog      'rm -vf *.log'
alias dscratch  'rm -vf CSCRATCH* BCSCRATCH* DSCRATCH*'
#
# Removes soft links from current directory
#
alias rmlinks         "find * -type l -maxdepth 0 -exec rm -vf {} ';'"
#
# Removes soft links from curent directory, and directories below it.
#
alias rm_all_links   "find * -type l -maxdepth 10 -exec rm -vf {} ';'"

