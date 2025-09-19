cd preprocessed

mkdir OCT_CytAssist_WT_713
cd OCT_CytAssist_WT_713

mkdir outs
cd outs

mkdir raw_feature_bc_matrix
cd raw_feature_bc_matrix


wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Fbarcodes%2Etsv%2Egz"
wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Ffeatures%2Etsv%2Egz"
wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Fmatrix%2Emtx%2Egz"

cd ..
mkdir spatial
cd spatial

wget -O tissue_lowres_image.png.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Ftissue%5Flowres%5Fimage%2Epng%2Egz"
gunzip tissue_lowres_image.png.gz

wget -O tissue_positions.csv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Ftissue%5Fpositions%2Ecsv%2Egz"
gunzip tissue_positions.csv.gz

wget -O scalefactors_json.json.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Fscalefactors%5Fjson%2Ejson%2Egz"
gunzip scalefactors_json.json.gz

wget -O detected_tissue_image.jpg.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8047nnn/GSM8047884/suppl/GSM8047884%5FOCT%5FCytAssist%5FWT%5F713%5Fdetected%5Ftissue%5Fimage%2Ejpg%2Egz"
gunzip detected_tissue_image.jpg.gz