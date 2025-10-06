cd preprocessed

mkdir FFPE_1_WT_460
cd FFPE_1_WT_460

mkdir outs
cd outs

mkdir raw_feature_bc_matrix
cd raw_feature_bc_matrix

wget -O barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Fmatrix%2Emtx%2Egz"

cd ..
mkdir spatial
cd spatial

wget -O tissue_lowres_image.png.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz

wget -O tissue_positions.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions.csv.gz

wget -O scalefactors_json.json.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz

wget -O detected_tissue_image.jpg.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8047868&format=file&file=GSM8047868%5FFFPE%5F1%5FWT%5F460%5Fdetected%5Ftissue%5Fimage%2Ejpg%2Egz"
gunzip detected_tissue_image.jpg.gz