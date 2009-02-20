#!/bin/bash

for dir in config/optflags/*/*; do
  if test "`ls ${dir} | grep all.standard`" = ""; then
    compilo=`dirname "${dir}"`
    compilo=`basename "${compilo}"`
    compilo=${compilo##*_}

    case ${compilo} in

      cc)
        echo 'CFLAGS_OPT="-O2"' > ${dir}/all.safe
        echo 'CFLAGS_OPT="-O2"' > ${dir}/all.standard
        echo 'CFLAGS_OPT="-O3"' > ${dir}/all.aggressive
        ;;

      cxx)
        echo 'CXXFLAGS_OPT="-O2"' > ${dir}/all.safe
        echo 'CXXFLAGS_OPT="-O2"' > ${dir}/all.standard
        echo 'CXXFLAGS_OPT="-O3"' > ${dir}/all.aggressive
        ;;

      fc)
        echo 'FCFLAGS_OPT="-O2"' > ${dir}/all.safe
        echo 'FCFLAGS_OPT="-O2"' > ${dir}/all.standard
        echo 'FCFLAGS_OPT="-O3"' > ${dir}/all.aggressive
        ;;

    esac
  fi
done

