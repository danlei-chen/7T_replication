set n=1 # change according to actual sample size

cd /cluster/iaslab/FSMAP/recon/

foreach ind (`seq 1 1 $n`)

#	set subj=`awk NR==$ind /cluster/iaslab/FSMAP/rsfmri/subj_g2.lst`
#	set run_no = `awk NR==$ind /cluster/iaslab/rsfmri/FSMAP/T1_scan_number.lst`
	set subj=`awk NR==$ind /cluster/iaslab/FSMAP/rsfmri/subj_g2_temp.lst`
	set run_no = `awk NR==$ind /cluster/iaslab/FSMAP/rsfmri/T1_scan_number_g2_temp.lst`

	mkdir $subj
	mkdir $subj/mri
	mkdir $subj/mri/orig

	mri_convert /cluster/iaslab/FSMAP/$subj/T1/$run_no/zreconT1EPI/miepi_MPRAGE_FOR_FS_r_bfc.nii.gz $subj/mri/orig/001.mgz

end
