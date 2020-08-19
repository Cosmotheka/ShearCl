#!/bin/bash

#nside=4096
#mem=5
#nc=24
nside=512
mem=0.5
nc=12
irot_0=0
nrot=1000
# Set to true if you want to run the cleanup stage
run_cleanup=true
# Set to true if you want to run the mapping stage
# (don't do this until cleanup has finished).
run_mapping=false
# Set to true if you want to run the C_ell stage
# (don't do this until mapping has finished).
run_cls_signal=false
# Set to true if you want to compute the bandpower windows
# (don't do this until the C_ell stage has finished).
run_windows=false
# Set to true if you want to compute the PSF null spectra
# and the spectra for rotated ellipticities.
# (don't do this until the mapping stage has finished).
run_cls_extra=false
# Set to true if you want to run the covariance stage
# (don't do this until the C_ell stage has finished).
run_cov_signal=false
# Set to true if you want to compute covariances for
# the PSF null tests (don't do this until the C_ell
# stage has finished).
run_cov_psf=false


mkdir -p /mnt/extraspace/damonge/S8z_data/outputs

echo "Catalog cleanup"
for b in {0..3}
do
    comment="binning_${b}"
    pyexec="addqueue -c ${comment} -q cmb -m 12 -n 1 /usr/bin/python3"
    comm="${pyexec} subsample.py --bin-number ${b}"
    echo ${comment}
    if [ $run_cleanup = true ] ; then
        ${comm}
    fi
done
echo " "
exit

echo "Mapping"
for b in {0..3}
do
    comment="mapping_${b}_ns${nside}_nr${nrot}"
    pyexec="addqueue -c ${comment} -q cmb -m 12 -n 1 /usr/bin/python3"
    comm="${pyexec} map.py --bin-number ${b} --nside ${nside} --nrot ${nrot}"
    echo ${comment}
    if [ $run_mapping = true ] ; then
        ${comm}
    fi
done
echo " "


echo "Cls - signal"
for b1 in {0..3}
do
    for b2 in {0..3}
    do
	if [[ ${b2} -lt ${b1} ]]; then
	    continue
	fi
	comment="ss_${b1}_${b2}_ns${nside}"
	pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
	comm="${pyexec} cls.py --bin-number ${b1} --bin-number-2 ${b2} --nside ${nside} --n-iter 0"
	echo ${comment}
        if [ $run_cls_signal = true ] ; then
            ${comm}
        fi
    done
done
echo " "


echo "Bandpower windows"
for b1 in {0..3}
do
    for b2 in {0..3}
    do
	if [[ ${b2} -lt ${b1} ]]; then
	    continue
	fi
	comment="win_${b1}_${b2}_ns${nside}"
	pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
	comm="${pyexec} windows.py --bin-number ${b1} --bin-number-2 ${b2} --nside ${nside}"
	echo ${comment}
        if [ $run_windows = true ] ; then
	    ${comm}
        fi
    done
done


echo "Cls - rotations and PSF"
for b in {0..3}
do
    #Rotations
    comment="rots_${b}_ns${nside}"
    pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
    comm="${pyexec} cls.py --bin-number ${b} --nside ${nside} --n-iter 0 --irot-0 ${irot_0} --irot-f ${nrot}"
    echo ${comment}
    if [ $run_cls_extra = true ] ; then
        ${comm}
    fi
    #PSF-x
    comment="psfX_${b}_ns${nside}"
    pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
    comm="${pyexec} cls.py --bin-number ${b} --nside ${nside} --n-iter 0 --is-psf-x"
    echo ${comment}
    if [ $run_cls_extra = true ] ; then
        ${comm}
    fi
    #PSF-a
    comment="psfA_${b}_ns${nside}"
    pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
    comm="${pyexec} cls.py --bin-number ${b} --nside ${nside} --n-iter 0 --is-psf-a"
    echo ${comment}
    if [ $run_cls_signal = true ] ; then
        ${comm}
    fi
done
echo " "


echo "Covs - signal"
icov=0
ia=0
for ba1 in {0..3}
do
    for ba2 in {0..3}
    do
	if [[ ${ba2} -lt ${ba1} ]]; then
	    continue
	fi

	ib=0
	for bb1 in {0..3}
	do
	    for bb2 in {0..3}
	    do
		if [[ ${bb2} -lt ${bb1} ]]; then
		    continue
		fi
		if [[ ${ib} -lt ${ia} ]]; then
		    ((ib++))
		    continue
		fi
		comment="cv_${ba1}${ba2}_${bb1}${bb2}_${ia}_${ib}_${icov}_ns${nside}"
		pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
		comm="${pyexec} covs.py --bin-a1 ${ba1} --bin-a2 ${ba2} --bin-b1 ${bb1} --bin-b2 ${bb2} --nside ${nside} --n-iter 0 --full-noise"

		echo ${comment}
                if [ $run_cov_signal = true ] ; then
                    ${comm}
                fi
		((ib++))
		((icov++))
	    done
	done
	((ia++))
    done
done
echo " "


echo "Covs - PSF"
for b in {0..3}
do
    comment="cv_xpsf_${b}_ns${nside}"
    pyexec="addqueue -c ${comment} -n 1x${nc} -s -q cmb -m ${mem} /usr/bin/python3"
    comm="${pyexec} covs_xPSF.py --bin-number ${b} --nside ${nside}"
    echo ${comment}
    if [ $run_cov_psf = true ] ; then
        ${comm}
    fi
done
echo " "

echo "Done!"
