#! /bin/bash
#grep -v '^#' psa64_pic_stripe_SEDfit_wFOM.txt | sort  -t' ' -k 8 -g  > psa64_pic_stripe_SEDfit_FOMsorted.txt
grep -v '^#' psa64_pic_stripe_perleycorr_SEDfits.txt | sort  -t' ' -k 10 -g  > \
psa64_pic_stripe_perleycorr_SEDfits_Pcorrsorted.txt

MYFILES=
i=0
for F in `cut -d' ' -f 1 psa64_pic_stripe_perleycorr_SEDfits_Pcorrsorted.txt`; 
do 
#print useful latex stuff
if [ $i -eq 15 ]
then
i=0
fi
if [ $i -eq 0 ]
then
echo "\begin{figure*}[htbp]"
echo "\begin{center}"
fi

MYFILES="${F}*MCMC.png $MYFILES"; 
MYFILE=`ls ${F}*MCMC.png`
echo "\includegraphics[width=2in]{plots/${MYFILE}}"   
((i++)) 

if [ $i -eq 15 ]
then
echo "\end{center}"
echo "\caption{fits of the next 16 sources, as described in Figure \ref{fig:src_spec1}."
echo "}\label{fig:src_spec1}"
echo "\end{figure*}"
echo "\clearpage"
i=0
fi
done
echo $MYFILES
open $MYFILES


