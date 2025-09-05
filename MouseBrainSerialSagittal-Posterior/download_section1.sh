cd preprocessed

mkdir section1
cd section1

mkdir outs
cd outs

wget https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior/V1_Mouse_Brain_Sagittal_Posterior_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior/V1_Mouse_Brain_Sagittal_Posterior_spatial.tar.gz

gunzip V1_Mouse_Brain_Sagittal_Posterior_raw_feature_bc_matrix.tar.gz
tar -xvf V1_Mouse_Brain_Sagittal_Posterior_raw_feature_bc_matrix.tar

gunzip V1_Mouse_Brain_Sagittal_Posterior_spatial.tar.gz 
tar -xvf V1_Mouse_Brain_Sagittal_Posterior_spatial.tar

# clean up
rm V1_Mouse_Brain_Sagittal_Posterior_raw_feature_bc_matrix.tar
rm V1_Mouse_Brain_Sagittal_Posterior_spatial.tar