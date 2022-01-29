home_dir=$(dirname $0)
cd $home_dir

if [ -d data_Zenodo ]; then
	echo cleaning up old data_Zenodo
	rm -rf ./data_Zenodo/*
else
	mkdir data_Zenodo
fi

wget https://zenodo.org/record/5919086/files/2022-01-29-data_Zenodo.tar.gz
tar -xvzf 2022-01-29-data_Zenodo.tar.gz -C data_Zenodo --strip-components=1
rm -f 2022-01-29-data_Zenodo.tar.gz