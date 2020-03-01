set n=1 # change according to actual sample size
set hpf=0.01
set lpf=0.08

cd /cluster/iaslab/FSMAP/

##################################################################################################################################################################
# CONFIGURATION SETUP (MANUALLY MODIFY AFTERWARDS)
#fcseed-config -vcsf -fcname vcsf.dat -fsd rest1 -pca -cfg rsfmri/vcsf1.config
#fcseed-config -vcsf -fcname vcsf.dat -fsd rest2 -pca -cfg rsfmri/vcsf2.config # THEN GO INTO FILE AND ADD 14 and 15 TO segidlists
#fcseed-config -wm -fcname wm.dat -fsd rest1 -pca -cfg rsfmri/wm1.config
#fcseed-config -wm -fcname wm.dat -fsd rest2 -pca -cfg rsfmri/wm2.config

#fcseed-config -vcsf -fcname vcsf_mean.dat -fsd rest1 -mean -cfg rsfmri/vcsf1_mean.config
#fcseed-config -wm -fcname wm_mean.dat -fsd rest1 -mean -cfg rsfmri/wm1_mean.config


# PREPROCESSING
foreach ind (`seq 1 1 $n`)

	set subj=`awk NR==$ind /cluster/iaslab/FSMAP/rsfmri/subj_g2_temp.lst` # change subj list
	set rest_dir='/cluster/iaslab/FSMAP/'$subj'/rest/'

	set first=`awk -v line=$ind -v col=1 'NR == '$ind' { print $1 }' < rsfmri/restingstate_scan_number_g2_temp.lst`
	set second=`awk -v line=$ind -v col=2 'NR == '$ind' { print $2 }' < rsfmri/restingstate_scan_number_g2_temp.lst`
	set third=`awk -v line=$ind -v col=3 'NR == '$ind' {print $3 }' < rsfmri/restingstate_scan_number_g2_temp.lst`
	set t1=`awk NR==$ind rsfmri/T1_scan_number_temp.lst`

	# SCAN NAME
	set fname=f_reorient

	# SLICE TIMING CORRECTION
	set TR=2.34 

	## SLICE TIMING CORRECTION
	foreach run ($first $second $third)

		#set run_dir=$rest_dir'/'$run

		## SPLIT THE TIMESERIES IN 3 SLABS
		fslroi  $run_dir/${fname}  $run_dir/${fname}_slab1 0 -1 0 -1 0 41 0 -1
		fslroi  $run_dir/${fname}  $run_dir/${fname}_slab2 0 -1 0 -1 41 41 0 -1
		fslroi  $run_dir/${fname}  $run_dir/${fname}_slab3 0 -1 0 -1 82 41 0 -1

		## PERFORM SLICE TIMING FOR EACH SLAB SEPARATELY
		slicetimer -i $run_dir/${fname}_slab1 -o $run_dir/${fname}_slab1_st --ocustom=/cluster/iaslab/FSMAP/scripts/sliceorder_123slsms3_asc_inter.txt  -v -r $TR
		slicetimer -i $run_dir/${fname}_slab2 -o $run_dir/${fname}_slab2_st --ocustom=/cluster/iaslab/FSMAP/scripts/sliceorder_123slsms3_asc_inter.txt  -v -r $TR
		slicetimer -i $run_dir/${fname}_slab3 -o $run_dir/${fname}_slab3_st --ocustom=/cluster/iaslab/FSMAP/scripts/sliceorder_123slsms3_asc_inter.txt  -v -r $TR

		## STICH 3 SLABS BACK TOGETHER
		fslmerge -z $run_dir/${fname}_st  $run_dir/${fname}_slab1_st $run_dir/${fname}_slab2_st $run_dir/${fname}_slab3_st
		rm $run_dir/*slab*

	end

	# MERGE 3 RUNS INTO 1 RUN
	fslmerge -t $subj/rest/all_merged $subj/rest/$first/f_reorient_st.nii.gz $subj/rest/$second/f_reorient_st.nii.gz $subj/rest/$third/f_reorient_st.nii.gz

	# COREGISTRATION TO T1EPI_MPRAGE (this step takes a while)
	epi_reg --epi=$subj/rest/all_merged --t1=$subj/T1/$t1/zreconT1EPI/miepi_MPRAGE_FOR_FS_r_bfc.nii.gz --t1brain=$subj/T1/$t1/zreconT1EPI/miepi_MPRAGE_FOR_FS_r_bfc.nii.gz --out=$subj/rest/all_merged_coreg.nii.gz

	# SPLIT COREGISTERED RUN TO 2 RUNS
	#mkdir $subj/rest1
	#mkdir $subj/rest1/001
	#mkdir $subj/rest2
	#mkdir $subj/rest2/001

	fslroi $subj/rest/all_merged_coreg.nii.gz $subj/rest1/001/f.nii.gz 0 384
	rm $subj/rest1/001/f.nii
	gunzip $subj/rest1/001/f.nii.gz
	fslroi $subj/rest/all_merged_coreg.nii.gz $subj/rest2/001/f_temp.nii.gz 384 384 

	#Compute the transformation that aligns the first timepoint of rest 2 to the first of rest1
	flirt -in $subj/rest2/001/f_temp.nii.gz -ref $subj/rest1/001/f.nii -out $subj/rest2/001/f1.nii.gz -omat $subj/rest2/001/f_1strun2_21strun1mat

	# Apply this transformation to the whole timecourse
	flirt -interp spline -in $subj/rest2/001/f_temp.nii.gz -ref $subj/rest1/001/f.nii -applyxfm -init $subj/rest2/001/f_1strun2_21strun1mat -out $subj/rest2/001/f.nii.gz 

	rm $subj/rest2/001/f_temp.nii.gz
	rm $subj/rest2/001/f1.nii.gz
	rm $subj/rest2/001/f.nii
	gunzip $subj/rest2/001/f.nii.gz

###########################

	# PREPROC-SESS AND NUISANCE SIGNAL EXTRACTION
	foreach set (1)

		# ACTUAL PREPROC-SESS LINE TO DO MOTION CORRECTION
		preproc-sess -s $subj -fsd rest$set -fwhm 0 -per-session -force
		ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/fmc.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.nii.gz
		ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/fmc.mcdat /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.mcdat
		ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/fmc.nii.gz.mclog /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.nii.gz.mclog
		ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/fmc.mat.aff12.1D /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.mat.aff12.1D

		# THIS IS DONE TO CREATE A PREPROC-SESS LOG TO TRICK FREESURFER INTO THINKING PREPROC IS DONE UP TO NORMALIZATION
		preproc-sess -s $subj -fsd rest$set -fwhm 0 -per-run -update -surface mni152.fnirt lhrh -mni305-1mm
		ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/fmc.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/f_st_mc.nii.gz
		mv /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni152.fnirt.lh.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni152.fnirt.lh_orig.nii.gz
		mv /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni152.fnirt.rh.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni152.fnirt.rh_orig.nii.gz
		mv /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni305.1mm.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni305.1mm_orig.nii.gz

		fcseed-sess -s $subj -cfg rsfmri/vcsf${set}_lat_ventricle.config -overwrite #add 5vcsf and aqueduct
		fcseed-sess -s $subj -cfg rsfmri/vcsf${set}_inf_lat_ventricle.config -overwrite
		fcseed-sess -s $subj -cfg rsfmri/vcsf${set}_choroid_plexus.config -overwrite
		fcseed-sess -s $subj -cfg rsfmri/vcsf${set}_3rd_ventricle.config -overwrite
		fcseed-sess -s $subj -cfg rsfmri/vcsf${set}_4th_ventricle.config -overwrite
		fcseed-sess -s $subj -cfg rsfmri/wm${set}_mean.config -overwrite
		fslmeants -i $subj/rest${set}/001/fmcpr.nii.gz -o $subj/rest${set}/001/aqueduct.dat -m $subj/aqueduct_mask_5vx.nii.gz
	end

	# CHECK STD IN COMBINED TIMESERIES	
		fslmerge -t $subj/combined_rest $subj/rest1/001/fmcpr.nii.gz $subj/rest2/001/fmcpr.nii.gz
		fslmaths $subj/combined_rest -Tstd $subj/combined_rest_STD
		rm $subj/combined_rest.nii.gz 

	# COMPUTE AQUEDUCT SIGNAL
		# fslview $subj/rest1/001/fmcpr.nii.gz $subj/combined_rest_STD -l Red-Yellow &

	## NUISANCE REGRESSION DONE IN MATLAB EXTERNAL TO THIS SCRIPT ##

	# BANDPASS FILTERING
	foreach set (1)
		# AFNI bandpass filtering
		# setenv PATH $PATH\:/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/AFNI/current/bin
		3dFourier -prefix $subj/rest$set/001/f_st_mc_regout_hp.nii.gz -highpass $hpf $subj/rest$set/001/f_st_mc_regout.nii.gz
		3dFourier -prefix $subj/rest$set/001/f_st_mc_regout_bp.nii.gz  -lowpass $lpf $subj/rest$set/001/f_st_mc_regout_hp.nii.gz
		rm $subj/rest$set/001/f_st_mc_regout_hp.nii.gz
	end

	## NORMALIZATION DONE EXTERNAL TO THIS SCRIPT ##

	foreach set (1)

		# CONVERSION TO FS ORIENTATION/DIMENSIONS 
		fslswapdim $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI.nii.gz x -z y $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA.nii.gz
		fslroi $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA.nii.gz $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop.nii.gz 14.5 151 14.5 151 16 186
		mri_convert $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop.nii.gz -oc -1.5 -16 9.5 $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc.nii.gz


		# DETREND
		fslmaths $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc.nii.gz -Tmean $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc_Tmean.nii.gz
		fslmaths $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc.nii.gz -sub $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc_Tmean.nii.gz $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_det.nii.gz
		ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_det.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni305.1mm.nii.gz

		rm $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA.nii.gz
		rm $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop.nii.gz
		rm $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc.nii.gz
		rm $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_LIA_crop_oc_Tmean.nii.gz

		# PROJECT TO LH/RH
		foreach hemi (lh rh)
			mri_vol2surf --src $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_det.nii.gz --srcreg identity --ref T1.mgz --regheader mni152.fnirt --hemi $hemi --o $subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_det_${hemi}.nii.gz
			ln -s /cluster/iaslab/FSMAP/$subj/rest$set/001/f_st_mc_regout_bp_m_2MNI_det_${hemi}.nii.gz /cluster/iaslab/FSMAP/$subj/rest$set/001/fmcpr.sm0.mni152.fnirt.${hemi}.nii.gz

		end
	end


end
