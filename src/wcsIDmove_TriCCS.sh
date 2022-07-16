#!/bin/bash
# Example 2022-04-11
# wcscheck_TriCCS_auto.sh cutedge wcs 5 18 5 18


DIR0="wcs"
DIR1="wcs_common"
mkdir -p ${DIR1}

mv hoge.txt ${DIR1}

while read x
do
  mv ${DIR0}/*${x}* ${DIR1}
done < idlist_0_common.txt

while read x
do
  mv ${DIR0}/*${x}* ${DIR1}
done < idlist_1_common.txt

while read x
do
  mv ${DIR0}/*${x}* ${DIR1}
done < idlist_2_common.txt
