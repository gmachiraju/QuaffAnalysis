#
# Makefile for Quaff Analysis on MinIon Oxford Nanopore Reads
# By Gautam Machiraju, Holmes Lab, UC Berkeley - Spring 2016
#

#-------------------EXAMPLE COMMANDS-------------------
# outputfile: inputfile.fa
#         quaff -some args inputfile.fa >outputfile

# outputs/%.txt: inputs/%.fa
#         quaff -some args inputs/$*.fa >outputs/$*.txt
#------------------------------------------------------

all: aligners rawdata d1 d2 d3 d4 d5 d1_train d2_train d3_train d4_train d5_train d1_align d2_align d3_align d4_align mA_train&align last_align bwa_align

clean:
	if [ -a quaff ]; then rm -rf quaff; fi;
	if [ -a marginAlign ]; then rm -rf marginAlign; fi;
	if [ -a last-737 ]; then rm -rf last-737; fi;
	if [ -a bwa-0.7.12 ]; then rm -rf bwa-0.7.12; fi;
	if [ -a __MACOSX ]; then rm -rf __MACOSX; fi;
	if [ -a D1 ]; then rm -rf D1; fi;
	if [ -a D2 ]; then rm -rf D2; fi;
	if [ -a D3 ]; then rm -rf D3; fi;

libraries:
	#- xcode tools (includes Clang, git), - Homebrew, - libz, libgsl

# Step 0:
# Dowload all the aligners 
aligners: 
	# Install and build quaff
	if [ -a quaff ]; then echo ""; else $ git clone https://github.com/ihh/quaff.git && cd quaff && make quaff && sudo make install; fi;

	# Install and build marginAlign
	if [ -a marginAlign ]; then echo ""; else $ git clone git://github.com/benedictpaten/marginAlign.git && cd marginAlign && $ git submodule update --init && make all; fi;

	# Install and build last
	if [ -a last-737 ]; then echo ""; else curl -o last-737.zip "http://last.cbrc.jp/last-737.zip" && unzip last-737.zip && cd last-737 && make && sudo make install && cd .. && rm last-737.zip; fi;

	# Install and build bwa
	# # curl -o bwa-0.7.12.tar.bz2 "http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2?r=&ts=1461874293&use_mirror=iweb"
	# # ^ will keep the zip file in initial folder  (* Nonfunctional!)
	if [ -a bwa-0.7.12 ]; then echo ""; else unzip bwa-0.7.12.zip; fi;


.PHONY: d4, d5

# Step 1-3:
# Download datafiles from EBI and place them into 4 main directories (D1-D4)

# Step 4:
# Split data files into 2 halves (training and testing) using file_splitter.py	
rawdata:
	mkdir rawdata

d1:
	mkdir D1

	# grabbing, unzipping
	mkdir D1/ERX708228
	curl -o ERX708228_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764952/ERR764952_1.fastq.gz"
	gunzip ERX708228_1.fastq.gz
	curl -o ERX708228_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764952/ERR764952_2.fastq.gz"
	gunzip ERX708228_2.fastq.gz
	mv *.fastq D1/ERX708228

	# merging, spliting
	$ python file_merger.py ERX708228.fastq D1/ERX708228
	mv ERX708228.fastq rawdata
	rm D1/ERX708228/ERX708228_1.fastq
	rm D1/ERX708228/ERX708228_2.fastq
	$ python file_splitter.py rawdata/ERX708228.fastq 2
	mv rawdata/*group.fastq D1/ERX708228

	#-------------------------------------------------------

	mkdir D1/ERX708229
	curl -o ERX708229_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764953/ERR764953_1.fastq.gz"
	gunzip ERX708229_1.fastq.gz
	curl -o ERX708229_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764953/ERR764953_2.fastq.gz"
	gunzip ERX708229_2.fastq.gz
	mv *.fastq D1/ERX708229

	$ python file_merger.py ERX708229.fastq D1/ERX708229
	mv ERX708229.fastq rawdata
	rm D1/ERX708229/ERX708229_1.fastq
	rm D1/ERX708229/ERX708229_2.fastq
	$ python file_splitter.py rawdata/ERX708229.fastq 2
	mv rawdata/*group.fastq D1/ERX708229

	#-------------------------------------------------------

	mkdir D1/ERX708230
	curl -o ERX708230_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764954/ERR764954_1.fastq.gz"
	gunzip ERX708230_1.fastq.gz
	curl -o ERX708230_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764954/ERR764954_2.fastq.gz"
	gunzip ERX708230_2.fastq.gz
	mv *.fastq D1/ERX708230

	$ python file_merger.py ERX708230.fastq D1/ERX708230
	mv ERX708230.fastq rawdata
	rm D1/ERX708230/ERX708230_1.fastq
	rm D1/ERX708230/ERX708230_2.fastq
	$ python file_splitter.py rawdata/ERX708230.fastq 2
	mv rawdata/*group.fastq D1/ERX708230

	#-------------------------------------------------------

	mkdir D1/ERX708231
	curl -o ERX708231_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764955/ERR764955_1.fastq.gz"
	gunzip ERX708231_1.fastq.gz
	curl -o ERX708231_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR764/ERR764955/ERR764955_2.fastq.gz"
	gunzip ERX708231_2.fastq.gz
	mv *.fastq D1/ERX708231

	$ python file_merger.py ERX708231.fastq D1/ERX708231
	mv ERX708231.fastq rawdata
	rm D1/ERX708231/ERX708231_1.fastq
	rm D1/ERX708231/ERX708231_2.fastq
	$ python file_splitter.py rawdata/ERX708231.fastq 2
	mv rawdata/*group.fastq D1/ERX708231


d2:
	mkdir D2

	mkdir D2/Forward
	mkdir D2/Reverse
	mkdir D2/HQ2D
	mkdir D2/Normal2D

	curl -o MG1655.fasta.tgz "ftp://climb.genomics.cn/pub/10.5524/100001_101000/100102/Ecoli_R7_CombinedFasta.tgz"
	tar -xvzf MG1655.fasta.tgz

	mv MG1655.fasta.tgz rawdata
	mv *.fasta rawdata
	
	$ python file_splitter.py rawdata/ForwardReads.fasta 2
	mv rawdata/*group.fasta D2/Forward
	$ python file_splitter.py rawdata/ReverseReads.fasta 2
	mv rawdata/*group.fasta D2/Reverse
	$ python file_splitter.py rawdata/HighQualityTwoDirectionReads.fasta 2
	mv rawdata/*group.fasta D2/HQ2D
	$ python file_splitter.py rawdata/NormalTwoDirectionReads.fasta 2
	mv rawdata/*group.fasta D2/Normal2D


d3:
	# finish in smiliar format to d1.
	mkdir D3
	
	# "MARC phase1b" (1-3)
	mkdir D3/ERX955609 # <----experiment accesssion number
	curl -o ERX955609_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968975/ERR968975.fastq.gz"
	gunzip ERX955609_1.fastq.gz
	curl -o ERX955609_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968975/ERR968975_1.fastq.gz"
	gunzip ERX955609_2.fastq.gz
	curl -o ERX955609_3.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968975/ERR968975_2.fastq.gz"
	gunzip ERX955609_3.fastq.gz
	mv *.fastq D3/ERX955609
	
	$ python file_merger.py ERX955609.fastq D3/ERX955609
	mv ERX955609.fastq rawdata
	rm D3/ERX955609/ERX955609_1.fastq
	rm D3/ERX955609/ERX955609_2.fastq
	rm D3/ERX955609/ERX955609_3.fastq
	$ python file_splitter.py rawdata/ERX955609.fastq 2
	mv rawdata/*group.fastq D3/ERX955609

	#-------------------------------------------------------
	# "MARC phase1b" (4-6)
	mkdir D3/ERX963580
	curl -o ERX963580_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968976/ERR968976.fastq.gz"
	gunzip ERX963580_1.fastq.gz
	curl -o ERX963580_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968976/ERR968976_1.fastq.gz"
	gunzip ERX963580_2.fastq.gz
	curl -o ERX963580_3.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968976/ERR968976_2.fastq.gz"
	gunzip ERX963580_3.fastq.gz
	mv *.fastq D3/ERX963580

	$ python file_merger.py ERX963580.fastq D3/ERX963580
	mv ERX963580.fastq rawdata
	rm D3/ERX963580/ERX963580_1.fastq
	rm D3/ERX963580/ERX963580_2.fastq
	rm D3/ERX963580/ERX963580_3.fastq
	$ python file_splitter.py rawdata/ERX963580.fastq 2
	mv rawdata/*group.fastq D3/ERX963580

	#-------------------------------------------------------
	# "Phase 1b Evaluation"
	mkdir D3/ERX978855
	curl -o ERX978855_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968972/ERR968972.fastq.gz"
	gunzip ERX978855_1.fastq.gz
	curl -o ERX978855_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968972/ERR968972_1.fastq.gz"
	gunzip ERX978855_2.fastq.gz
	curl -o ERX978855_3.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR968/ERR968972/ERR968972_2.fastq.gz"
	gunzip ERX978855_3.fastq.gz
	mv *.fastq D3/ERX978855

	$ python file_merger.py ERX978855.fastq D3/ERX978855
	mv ERX978855.fastq rawdata
	rm D3/ERX978855/ERX978855_1.fastq
	rm D3/ERX978855/ERX978855_2.fastq
	rm D3/ERX978855/ERX978855_3.fastq
	$ python file_splitter.py rawdata/ERX978855.fastq 2
	mv rawdata/*group.fastq D3/ERX978855

d4:
	# D4 will be given in makefile in folder D4 package to speed up download.
	# Need to move it to raw data and split it up, then move split files to D4...
	unzip D4/MAP006-1.fasta.zip
	mv D4/MAP006-1.fasta.zip rawdata/MAP006-1.fasta.zip
	mv MAP006-1.fasta rawdata
	$ python file_splitter.py rawdata/MAP006-1.fasta 2
	mv rawdata/*group.fasta D4
	

d5:
	# Already included in the package as well, as D5.
	# Nothing to be done.
	unzip D5/U00096.3.fasta.zip
	mv D5/U00096.3.fasta.zip rawdata/U00096.3.fasta.zip
	mv U00096.3.fasta D5


# Step 5:
# Run "quaff train" to fit the model to the E.coli genome (D5), and save the trained 
# parameter JSON files for each dataset's training half (D1-D4)
# 		Do this first on 1% of the read data (using -maxreadmb), then on 100%. Repeat this 
# 		for order-1 and order-2 models (using -order). So that is a total of 24 experiments: 
# 		(datasets D1 through D4) * (1% and 100%) * (order-0, order-1 and order-2) = 4*2*3 = 24
d1_train:
	#renaming
	cp D1/ERX708228/ERX708228_1_group.fastq D1/ERX708228/ERX708228_1.fastq
	cp D1/ERX708228/ERX708228_2_group.fastq D1/ERX708228/ERX708228_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D1/ERX708228/ERX708228_1.fastq 100

	# D1-228, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708228/ERX708228_1_1_group.fastq -order 0  >paramsD1-228-1-0.json
	# D1-228, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708228/ERX708228_1.fastq -order 0  >paramsD1-228-100-0.json

	# D1-228, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708228/ERX708228_1_1_group.fastq -order 1  >paramsD1-228-1-1.json
	# D1-228, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708228/ERX708228_1.fastq -order 1  >paramsD1-228-100-1.json

	# D1-228, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D1/ERX708228/ERX708228_1_1_group.fastq -order 2  >paramsD1-228-1-2.json
	# D1-228, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D1/ERX708228ERX708228_1.fastq -order 2  >paramsD1-228-100-2.json

	#-----------------------------------------------------------------------------------

	#renaming
	cp D1/ERX708229/ERX708229_1_group.fastq D1/ERX708229/ERX708229_1.fastq
	cp D1/ERX708229/ERX708229_2_group.fastq D1/ERX708229/ERX708229_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D1/ERX708229/ERX708229_1.fastq 100

	# D1-229, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708229/ERX708229_1_1_group.fastq -order 0  >paramsD1-229-1-0.json
	# D1-229, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708229/ERX708229_1.fastq -order 0  >paramsD1-229-100-0.json

	# D1-229, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708229/ERX708229_1_1_group.fastq -order 1  >paramsD1-229-1-1.json
	# D1-229, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708229/ERX708229_1.fastq -order 1  >paramsD1-229-100-1.json

	# D1-229, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D1/ERX708229/ERX708229_1_1_group.fastq -order 2  >paramsD1-229-1-2.json
	# D1-229, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D1/ERX708229/ERX708229_1.fastq -order 2  >paramsD1-229-100-2.json


	#-----------------------------------------------------------------------------------

	#renaming
	cp D1/ERX708230/ERX708230_1_group.fastq D1/ERX708230/ERX708230_1.fastq
	cp D1/ERX708230/ERX708230_2_group.fastq D1/ERX708230/ERX708230_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D1/ERX708230/ERX708230_1.fastq 100

	# D1-230, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708230/ERX708230_1_1_group.fastq -order 0  >paramsD1-230-1-0.json
	# D1-230, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708230/ERX708230_1.fastq -order 0  >paramsD1-230-100-0.json

	# D1-230, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708230/ERX708230_1_1_group.fastq -order 1  >paramsD1-230-1-1.json
	# D1-230, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708230/ERX708230_1.fastq -order 1  >paramsD1-230-100-1.json

	# D1-230, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D1/ERX708230/ERX708230_1_1_group.fastq -order 2  >paramsD1-230-1-2.json
	# D1-230, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D1/ERX708230/ERX708230_1.fastq -order 2  >paramsD1-230-100-2.json


	#-----------------------------------------------------------------------------------

	#renaming
	cp D1/ERX708231/ERX708231_1_group.fastq D1/ERX708231/ERX708231_1.fastq
	cp D1/ERX708231/ERX708231_2_group.fastq D1/ERX708231/ERX708231_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D1/ERX708231/ERX708231_1.fastq 100

	# D1-231, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708231/ERX708231_1_1_group.fastq -order 0 >paramsD1-231-1-0.json
	# D1-231, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D1/ERX708231/ERX708231_1.fastq -order 0 >paramsD1-231-100-0.json

	# D1-231, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708231/ERX708231_1_1_group.fastq -order 1 >paramsD1-231-1-1.json
	# D1-231, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D1/ERX708231/ERX708231_1.fastq -order 1 >paramsD1-231-100-1.json

	# D1-231, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D1/ERX708231/ERX708231_1_1_group.fastq -order 2 >paramsD1-231-1-2.json
	# D1-231, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D1/ERX708231/ERX708231_1.fastq -order 2 >paramsD1-231-100-2.json


d2_train:
	cp D2/Forward/ForwardReads_1_group.fasta D2/Forward/ForwardReads_1.fasta
	cp D2/Forward/ForwardReads_2_group.fasta D2/Forward/ForwardReads_2.fasta

	$ python file_splitter.py D2/Forward/ForwardReads_1.fasta 100

	# D2-For, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/Forward/ForwardReads_1_1_group.fasta -order 0 >paramsD2-For-1-0.json
	# D2-For, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/Forward/ForwardReads_1.fasta -order 0 >paramsD2-For-100-0.json

	# D2-For, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D2/Forward/ForwardReads_1_1_group.fasta -order 1 >paramsD2-For-1-1.json
	# D2-For, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D2/Forward/ForwardReads_1.fasta -order 1 >paramsD2-For-100-1.json

	# D2-For, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D2/Forward/ForwardReads_1_1_group.fasta -order 2 >paramsD2-For-1-2.json
	# D2-For, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D2/Forward/ForwardReads_1.fasta -order 2 >paramsD2-For-100-2.json

	#-----------------------------------------------------------------------------------

	cp D2/HQ2D/HighQualityTwoDirectionReads_1_group.fasta D2/HQ2D/HQ2DReads_1.fasta
	cp D2/HQ2D/HighQualityTwoDirectionReads_2_group.fasta D2/HQ2D/HQ2DReads_2.fasta

	$ python file_splitter.py D2/HQ2D/HQ2DReads_1.fasta 100

	# D2-HQ2D, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/HQ2D/HQ2DReads_1_1_group.fasta -order 0 >paramsD2-HQ2D-1-0.json
	# D2-HQ2D, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/HQ2D/HQ2DReads_1.fasta -order 0 >paramsD2-HQ2D-100-0.json

	# D2-HQ2D, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D2/HQ2D/HQ2DReads_1_1_group.fasta -order 1 >paramsD2-HQ2D-1-1.json
	# D2-HQ2D, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D2/HQ2D/HQ2DReads_1.fasta -order 1 >paramsD2-HQ2D-100-1.json

	# D2-HQ2D, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D2/HQ2D/HQ2DReads_1_1_group.fasta -order 2 >paramsD2-HQ2D-1-2.json
	# D2-HQ2D, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D2/HQ2D/HQ2DReads_1.fasta -order 2 >paramsD2-HQ2D-100-2.json

	#-----------------------------------------------------------------------------------

	cp D2/Normal2D/NormalTwoDirectionReads_1_group.fasta D2/Normal2D/Nr2DReads_1.fasta
	cp D2/Normal2D/NormalTwoDirectionReads_2_group.fasta D2/Normal2D/Nr2DReads_2.fasta

	$ python file_splitter.py D2/Normal2D/Nr2DReads_1.fasta 100

	# D2-Nr2D, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/Normal2D/Nr2DReads_1_1_group.fasta -order 0 >paramsD2-Nr2D-1-0.json
	# D2-Nr2D, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/Normal2D/Nr2DReads_1.fasta -order 0 >paramsD2-Nr2D-100-0.json

	# D2-Nr2D, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D2/Normal2D/Nr2DReads_1_1_group.fasta -order 1 >paramsD2-Nr2D-1-1.json
	# D2-Nr2D, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D2/Normal2D/Nr2DReads_1.fasta -order 1 >paramsD2-Nr2D-100-1.json

	# D2-Nr2D, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D2/Normal2D/Nr2DReads_1_1_group.fasta -order 2 >paramsD2-Nr2D-1-2.json
	# D2-Nr2D, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D2/Normal2D/Nr2DReads_1.fasta -order 2 >paramsD2-Nr2D-100-2.json

	#-----------------------------------------------------------------------------------

	cp D2/Reverse/ReverseReads_1_group.fasta D2/Reverse/ReverseReads_1.fasta
	cp D2/Reverse/ReverseReads_2_group.fasta D2/Reverse/ReverseReads_2.fasta

	$ python file_splitter.py D2/Reverse/ReverseReads_1.fasta 100

	# D2-Rev, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/Reverse/ReverseReads_1_1_group.fasta -order 0 >paramsD2-Rev-1-0.json
	# D2-Rev, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D2/Reverse/ReverseReads_1.fasta -order 0 >paramsD2-Rev-100-0.json

	# D2-Rev, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D2/Reverse/ReverseReads_1_1_group.fasta -order 1 >paramsD2-Rev-1-1.json
	# D2-Rev, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D2/Reverse/ReverseReads_1.fasta -order 1 >paramsD2-Rev-100-1.json

	# D2-Rev, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D2/Reverse/ReverseReads_1_1_group.fasta -order 2 >paramsD2-Rev-1-2.json
	# D2-Rev, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D2/Reverse/ReverseReads_1.fasta -order 2 >paramsD2-Rev-100-2.json


d3_train:
	#renaming
	cp D3/ERX955609/ERX955609_1_group.fastq D3/ERX955609/ERX955609_1.fastq
	cp D3/ERX955609/ERX955609_2_group.fastq D3/ERX955609/ERX955609_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D3/ERX955609/ERX955609_1.fastq 100

	# D3-609, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D3/ERX955609/ERX955609_1_1_group.fastq -order 0 >paramsD3-609-1-0.json
	# D3-609, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D3/ERX955609/ERX955609_1.fastq -order 0 >paramsD3-609-100-0.json

	# D3-609, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D3/ERX955609/ERX955609_1_1_group.fastq -order 1 >paramsD3-609-1-1.json
	# D3-609, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D3/ERX955609/ERX955609_1.fastq -order 1 >paramsD3-609-100-1.json

	# D3-609, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D3/ERX955609/ERX955609_1_1_group.fastq -order 2 >paramsD3-609-1-2.json
	# D3-609, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D3/ERX955609/ERX955609_1.fastq -order 2 >paramsD3-609-100-2.json

	#-----------------------------------------------------------------------------------

	#renaming
	cp D3/ERX963580/ERX963580_1_group.fastq D3/ERX963580/ERX963580_1.fastq
	cp D3/ERX963580/ERX963580_2_group.fastq D3/ERX963580/ERX963580_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D3/ERX963580/ERX963580_1.fastq 100

	# D3-580, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D3/ERX963580/ERX963580_1_1_group.fastq -order 0 >paramsD3-580-1-0.json
	# D3-580, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D3/ERX963580/ERX963580_1.fastq -order 0 >paramsD3-580-100-0.json

	# D3-580, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D3/ERX963580/ERX963580_1_1_group.fastq -order 1 >paramsD3-580-1-1.json
	# D3-580, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D3/ERX963580/ERX963580_1.fastq -order 1 >paramsD3-580-100-1.json

	# D3-580, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D3/ERX963580/ERX963580_1_1_group.fastq -order 2 >paramsD3-580-1-2.json
	# D3-580, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D3/ERX963580/ERX963580_1.fastq -order 2 >paramsD3-580-100-2.json

	#-----------------------------------------------------------------------------------

	#renaming
	cp D3/ERX978855/ERX978855_1_group.fastq D3/ERX978855/ERX978855_1.fastq
	cp D3/ERX978855/ERX978855_2_group.fastq D3/ERX978855/ERX978855_2.fastq

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D3/ERX978855/ERX978855_1.fastq 100

	# D3-855, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D3/ERX978855/ERX978855_1_1_group.fastq -order 0 >paramsD3-855-1-0.json
	# D3-855, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D3/ERX978855/ERX978855_1.fastq -order 0 >paramsD3-855-100-0.json

	# D3-855, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D3/ERX978855/ERX978855_1_1_group.fastq -order 1 >paramsD3-855-1-1.json
	# D3-855, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D3/ERX978855/ERX978855_1.fastq -order 1 >paramsD3-855-100-1.json

	# D3-855, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D3/ERX978855/ERX978855_1_1_group.fastq -order 2 >paramsD3-855-1-2.json
	# D3-855, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D3/ERX978855/ERX978855_1.fastq -order 2 >paramsD3-855-100-2.json


d4_train: 
	#renaming
	cp D4/MAP006-1_1_group.fasta D4/MAP006-1_1.fasta
	cp D4/MAP006-1_2_group.fasta D4/MAP006-1_2.fasta

	#used this instead of maxreadmb -- need percentage..
	$ python file_splitter.py D4/MAP006-1_1.fasta 100

	# D4-MAP, 1%, ORDER 0 
	quaff train D5/U00096.3.fasta D4/MAP006-1_1_1_group.fasta -order 0 >paramsD4-MAP-1-0.json
	# D4-MAP, 100%, ORDER 0 
	quaff train D5/U00096.3.fasta D4/MAP006-1_1.fasta -order 0 >paramsD4-MAP-100-0.json

	# D4-MAP, 1%, ORDER 1
	quaff train D5/U00096.3.fasta D4/MAP006-1_1_1_group.fasta -order 1 >paramsD4-MAP-1-1.json
	# D4-MAP, 100%, ORDER 1
	quaff train D5/U00096.3.fasta D4/MAP006-1_1.fasta -order 1 >paramsD4-MAP-100-1.json

	# D4-MAP, 1%, ORDER 2 
	quaff train D5/U00096.3.fasta D4/MAP006-1_1_1_group.fasta -order 2 >paramsD4-MAP-1-2.json
	# D4-MAP, 100%, ORDER 2
	quaff train D5/U00096.3.fasta D4/MAP006-1_1.fasta -order 2 >paramsD4-MAP-100-2.json


# Step 6:
# Use "quaff align" to map the corresponding test sets (the other half of the dataset) back to 
# the E.coli genome for each of the trained parameter files. 
# 		Collect enough data to figure out what proportion of test-set reads align to the genome. 
#		(Collecting all the alignments should be enough.)
d1_align:
	# D1-228, 1%, ORDER 0
	quaff align D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq -params paramsD1-228-1-0.json >alignD1-228-1-0.stockholm 
	# D1-228, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq -params paramsD1-228-100-0.json >alignD1-228-100-0.stockholm
	# D1-228, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq -params paramsD1-228-1-1.json >alignD1-228-1-1.stockholm
	# D1-228, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq -params paramsD1-228-100-1.json >alignD1-228-100-1.stockholm
	# D1-228, 1%, ORDER 2
	quaff align D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq -params paramsD1-228-1-2.json >alignD1-228-1-2.stockholm 
	# D1-228, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq -params paramsD1-228-100-2.json >alignD1-228-100-2.stockholm

	# D1-229, 1%, ORDER 0
	quaff align D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq -params paramsD1-229-1-0.json >alignD1-229-1-0.stockholm
	# D1-229, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq -params paramsD1-229-100-0.json >alignD1-229-100-0.stockholm
	# D1-229, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq -params paramsD1-229-1-1.json >alignD1-229-1-1.stockholm
	# D1-229, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq -params paramsD1-229-100-1.json >alignD1-229-100-1.stockholm
	# D1-229, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq -params paramsD1-229-1-2.json >alignD1-229-1-2.stockholm
	# D1-229, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq -params paramsD1-229-100-2.json >alignD1-229-100-2.stockholm
 
	# D1-230, 1%, ORDER 0
	quaff align D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq -params paramsD1-230-1-0.json >alignD1-230-1-0.stockholm 
	# D1-230, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq -params paramsD1-230-100-0.json >alignD1-230=100-0.stockholm
	# D1-230, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq -params paramsD1-230-1-1.json >alignD1-230-1-1.stockholm
	# D1-230, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq -params paramsD1-230-100-1.json >alignD1-230-100-1.stockholm
	# D1-230, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq -params paramsD1-230-1-2.json >alignD1-230-1-2.stockholm
	# D1-230, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq -params paramsD1-230-100-2.json >alignD1-230-100-2.stockholm

	# D1-231, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq -params paramsD1-231-1-0.json >alignD1-231-1-0.stockholm
	# D1-231, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq -params paramsD1-231-100-0.json >alignD1-231-100-0.stockholm
	# D1-231, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq -params paramsD1-231-1-1.json >alignD1-231-1-1.stockholm
	# D1-231, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq -params paramsD1-231-100-1.json >alignD1-231-100-1.stockholm
	# D1-231, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq -params paramsD1-231-1-2.json >alignD1-231-1-2.stockholm
	# D1-231, 100%, ORDER 2
	quaff align D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq -params paramsD1-231-100-2.json >alignD1-231-100-2.stockholm 


d2_align:
	# D2-For, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta -params paramsD2-For-1-0.json >alignD2-For-1-0.stockholm
	# D2-For, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta -params paramsD2-For-100-0.json >alignD2-For-100-0.stockholm
	# D2-For, 1%, ORDER 1 
	quaff align D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta -params paramsD2-For-1-1.json >alignD2-For-1-1.stockholm
	# D2-For, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta -params paramsD2-For-100-1.json >alignD2-For-100-1.stockholm
	# D2-For, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta -params paramsD2-For-1-2.json >alignD2-For-1-2.stockholm
	# D2-For, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta -params paramsD2-For-100-2.json >alignD2-For-100-2.stockholm

	# D2-HQ2D, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta -params paramsD2-HQ2D-1-0.json >alignD2-HQ2D=1-0.stockholm
	# D2-HQ2D, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta -params paramsD2-HQ2D-100-0.json >alignD2-HQ2D-100-0.stockholm
	# D2-HQ2D, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta -params paramsD2-HQ2D-1-1.json >alignD2-HQ2D-1-1.stockholm
	# D2-HQ2D, 100%, ORDER 1
	quaff align D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta -params paramsD2-HQ2D-100-1.json >alignD2-HQ2D-100-1.stockholm
	# D2-HQ2D, 1%, ORDER 2
	quaff align D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta -params paramsD2-HQ2D-1-2.json >alignD2-HQ2D-1-2.stockholm
	# D2-HQ2D, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta -params paramsD2-HQ2D-100-2.json >alignD2-HQ2D-100-2.stockholm

	# D2-Nr2D, 1%, ORDER 0
	quaff align D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta -params paramsD2-Nr2D-1-0.json >alignD2-Nr2D-1-0.stockholm
	# D2-Nr2D, 100%, ORDER 0
	quaff align D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta -params paramsD2-Nr2D-100-0.json >alignD2-Nr2D-100-0.stockholm
	# D2-Nr2D, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta -params paramsD2-Nr2D-1-1.json >alignD2-Nr2D-1-1.stockholm
	# D2-Nr2D, 100%, ORDER 1
	quaff align D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta -params paramsD2-Nr2D-100-1.json >alignD2-Nr2D-100-1.stockholm
	# D2-Nr2D, 1%, ORDER 2
	quaff align D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta -params paramsD2-Nr2D-1-2.json >alignD2-Nr2D-1-2.stockholm
	# D2-Nr2D, 100%, ORDER 2
	quaff align D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta -params paramsD2-Nr2D-100-2.json >alignD2-Nr2D-100-2.stockholm

	# D2-Rev, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta -params paramsD2-Rev-1-0.json >alignD2-Rev-1-0.stockholm
	# D2-Rev, 100%, ORDER 0
	quaff align D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta -params paramsD2-Rev-100-0.json >alignD2-Rev-100-0.stockholm
	# D2-Rev, 1%, ORDER 1 
	quaff align D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta -params paramsD2-Rev-1-1.json >alignD2-Rev-1-1.stockholm
	# D2-Rev, 100%, ORDER 1
	quaff align D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta -params paramsD2-Rev-100-1.json >alignD2-Rev-100-1.stockholm
	# D2-Rev, 1%, ORDER 2
	quaff align D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta -params paramsD2-Rev-1-2.json >alignD2-Rev-1-2.stockholm
	# D2-Rev, 100%, ORDER 2
	quaff align D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta -params paramsD2-Rev-100-2.json >alignD2-Rev-100-2.stockholm

d3_align:
	# D3-609, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq -params paramsD3-609-1-0.json >alignD3-609-1-0.stockholm
	# D3-609, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq -params paramsD3-609-100-0.json >alignD3-609-100-0.stockholm
	# D3-609, 1%, ORDER 1 
	quaff align D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq -params paramsD3-609-1-1.json >alignD3-609-1-1.stockholm
	# D3-609, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq -params paramsD3-609-100-1.json >alignD3-609-100-1.stockholm
	# D3-609, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq -params paramsD3-609-1-2.json >alignD3-609-1-2.stockholm
	# D3-609, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq -params paramsD3-609-100-2.json >alignD3-609-100-2.stockholm

	# D3-580, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq -params paramsD3-580-1-0.json >alignD3-580-1-0.stockholm
	# D3-580, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq -params paramsD3-580-100-0.json >alignD3-580-100-0.stockholm
	# D3-580, 1%, ORDER 1 
	quaff align D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq -params paramsD3-580-1-1.json >alignD3-580-1-1.stockholm
	# D3-580, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq -params paramsD3-580-100-1.json >alignD3-580-100-1.stockholm
	# D3-580, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq -params paramsD3-580-1-2.json >alignD3-580-1-2.stockholm
	# D3-580, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq -params paramsD3-580-100-2.json >alignD3-580-100-2.stockholm

	# D3-855, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq -params paramsD3-855-1-0.json >alignD3-855-1-0.stockholm
	# D3-855, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq -params paramsD3-855-100-0.json >alignD3-855-100-0.stockholm
	# D3-855, 1%, ORDER 1 
	quaff align D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq -params paramsD3-855-1-1.json >alignD3-855-1-1.stockholm
	# D3-855, 100%, ORDER 1 
	quaff align D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq -params paramsD3-855-100-1.json >alignD3-855-100-1.stockholm
	# D3-855, 1%, ORDER 2 
	quaff align D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq -params paramsD3-855-1-2.json >alignD3-855-1-2.stockholm
	# D3-855, 100%, ORDER 2 
	quaff align D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq -params paramsD3-855-100-2.json >alignD3-855-100-2.stockholm


d4_align:
	# D4-MAP, 1%, ORDER 0 
	quaff align D5/U00096.3.fasta D4/MAP006-1_2.fasta -params paramsD4-MAP-1-0.json >alignD4-MAP.stockholm
	# D4-MAP, 100%, ORDER 0 
	quaff align D5/U00096.3.fasta D4/MAP006-1_2.fasta -params paramsD4-MAP-100-0.json >alignD4-MAP.stockholm
	# D4-MAP, 1%, ORDER 1
	quaff align D5/U00096.3.fasta D4/MAP006-1_2.fasta -params paramsD4-MAP-1-1.json >alignD4-MAP.stockholm
	# D4-MAP, 100%, ORDER 1
	quaff align D5/U00096.3.fasta D4/MAP006-1_2.fasta -params paramsD4-MAP-100-1.json >alignD4-MAP.stockholm
	# D4-MAP, 1%, ORDER 2
	quaff align D5/U00096.3.fasta D4/MAP006-1_2.fasta -params paramsD4-MAP-1-2.json >alignD4-MAP.stockholm
	# D4-MAP, 100%, ORDER 2
	quaff align D5/U00096.3.fasta D4/MAP006-1_2.fasta -params paramsD4-MAP-100-2.json >alignD4-MAP.stockholm


#--------------OTHER ALIGNER DATA GENERATION--------------

# Step 7:
# Repeat steps 5 and 6 for marginAlign (fitting it to each of datasets D1-D4 at both 1% and 100% of the full dataset)
mA_train_align:
	# Syntax:
	# reg-align:  marginAlign input.fastq reference.fasta output.sam --jobTree ./jobTree
	# train:      marginAlign input.fastq reference.fasta output.sam --em --outputModel output.hmm --jobTree ./jobTree
	# align:      marginAlign input.fastq reference.fasta output.sam --inputModel input.hmm --jobTree ./jobTree
	
	# D1 train: 
	#-----------
	# D1-228, 1%
	marginAlign/marginAlign D1/ERX708228/ERX708228_1_1_group.fastq D5/U00096.3.fasta paramsD1-228-1.sam --em --outputModel modelD1-228-1.hmm --jobTree ./jobTree228-1
	# D1-228, 100%
	marginAlign/marginAlign D1/ERX708228/ERX708228_1.fastq D5/U00096.3.fasta paramsD1-228-100.sam --em --outputModel modelD1-228-100.hmm --jobTree ./jobTree228-100

	# D1-229, 1%
	marginAlign/marginAlign D1/ERX708229/ERX708229_1_1_group.fastq D5/U00096.3.fasta paramsD1-229-1.sam --em --outputModel modelD1-229-1.hmm --jobTree ./jobTree229-1
	# D1-229, 100%
	marginAlign/marginAlign D1/ERX708229/ERX708229_1.fastq D5/U00096.3.fasta paramsD1-229-100.sam --em --outputModel modelD1-229-100.hmm --jobTree ./jobTree229-100

	# D1-230, 1%
	marginAlign/marginAlign D1/ERX708230/ERX708230_1_1_group.fastq D5/U00096.3.fasta paramsD1-230-1.sam --em --outputModel modelD1-230-1.hmm --jobTree ./jobTree230-1
	# D1-230, 100%
	marginAlign/marginAlign D1/ERX708230/ERX708230_1.fastq D5/U00096.3.fasta paramsD1-230-100.sam --em --outputModel modelD1-230-100.hmm --jobTree ./jobTree230-100

	# D1-231, 1%
	marginAlign/marginAlign D1/ERX708231/ERX708231_1_1_group.fastq D5/U00096.3.fasta paramsD1-231-1.sam --em --outputModel modelD1-231-1.hmm --jobTree ./jobTree231-1
	# D1-231, 100%
	marginAlign/marginAlign D1/ERX708231/ERX708231_1.fastq D5/U00096.3.fasta paramsD1-231-100.sam --em --outputModel modelD1-231-100.hmm --jobTree ./jobTree231-100
    

	# D1 align: 
	#-----------     
	marginAlign/marginAlign D1/ERX708228/ERX708228_2.fastq D5/U00096.3.fasta outD1-228-1.sam --inputModel modelD1-228-1.hmm --jobTree ./jobTree228-1
	marginAlign/marginAlign D1/ERX708228/ERX708228_2.fastq D5/U00096.3.fasta outD1-228-1.sam --inputModel modelD1-228-100.hmm --jobTree ./jobTree228-100

	marginAlign/marginAlign D1/ERX708229/ERX708229_2.fastq D5/U00096.3.fasta outD1-229-100.sam --inputModel modelD1-229-1.hmm --jobTree ./jobTree229-1
	marginAlign/marginAlign D1/ERX708229/ERX708229_2.fastq D5/U00096.3.fasta outD1-229-100.sam --inputModel modelD1-229-100.hmm --jobTree ./jobTree229-100

	marginAlign/marginAlign D1/ERX708230/ERX708230_2.fastq D5/U00096.3.fasta outD1-230-100.sam --inputModel modelD1-230-1.hmm --jobTree ./jobTree230-1
	marginAlign/marginAlign D1/ERX708230/ERX708230_2.fastq D5/U00096.3.fasta outD1-230-100.sam --inputModel modelD1-230-100.hmm --jobTree ./jobTree230-100

	marginAlign/marginAlign D1/ERX708231/ERX708231_2.fastq D5/U00096.3.fasta outD1-231-100.sam --inputModel modelD1-231-1.hmm --jobTree ./jobTree231-1
	marginAlign/marginAlign D1/ERX708231/ERX708231_2.fastq D5/U00096.3.fasta outD1-231-100.sam --inputModel modelD1-231-100.hmm --jobTree ./jobTree231-100


	# D2 train: 
	#-----------
	# D2-For, 1%
	marginAlign/marginAlign D2/Forward/ForwardReads_1_1_group.fasta D5/U00096.3.fasta paramsD2-For-1.sam --em --outputModel modelD2-For-1.hmm --jobTree ./jobTree-f1
	# D2-For, 100%
	marginAlign/marginAlign D2/Forward/ForwardReads_1.fasta D5/U00096.3.fasta paramsD2-For-100.sam --em --outputModel modelD2-For-100.hmm --jobTree ./jobTree-f100

	# D2-HQ2D, 1%
	marginAlign/marginAlign D2/HQ2D/HQ2DReads_1_1_group.fasta D5/U00096.3.fasta paramsD2-HQ2D-1.sam --em --outputModel modelD2-HQ2D-1.hmm --jobTree ./jobTree-h1
	# D2-HQ2D, 100%
	marginAlign/marginAlign D2/HQ2D/HQ2DReads_1.fasta D5/U00096.3.fasta paramsD2-HQ2D-100.sam --em --outputModel modelD2-HQ2D-100.hmm --jobTree ./jobTree-h100

	# D2-Nr2D, 1%
	marginAlign/marginAlign D2/Normal2D/Nr2DReads_1_1_group.fasta D5/U00096.3.fasta paramsD2-Nr2D-1.sam --em --outputModel modelD2-Nr2D-1.hmm --jobTree ./jobTree-n1
	# D2-Nr2D, 100%
	marginAlign/marginAlign D2/Normal2D/Nr2DReads_1.fasta D5/U00096.3.fasta paramsD2-Nr2D-100.sam --em --outputModel modelD2-Nr2D-100.hmm --jobTree ./jobTree-n100
	
	# D2-Rev, 1%	
	marginAlign/marginAlign D2/Reverse/ReverseReads_1_1_group.fasta D5/U00096.3.fasta paramsD2-Rev-1.sam --em --outputModel modelD2-Rev-1.hmm --jobTree ./jobTree-r1
	# D2-Rev, 100%	
	marginAlign/marginAlign D2/Reverse/ReverseReads_1.fasta D5/U00096.3.fasta paramsD2-Rev-100.sam --em --outputModel modelD2-Rev-100.hmm --jobTree ./jobTree-r100


	# D2 align: 
	#-----------     
	marginAlign/marginAlign D2/Forward/ForwardReads_2.fasta D5/U00096.3.fasta outD2-For-1.sam --inputModel modelD2-For-1.hmm --jobTree ./jobTree-f1
	marginAlign/marginAlign D2/Forward/ForwardReads_2.fasta D5/U00096.3.fasta outD2-For-100.sam --inputModel modelD2-For-100.hmm --jobTree ./jobTree-f100

	marginAlign/marginAlign D2/HQ2D/HQ2DReads_2.fasta D5/U00096.3.fasta outD2-HQ2D-1.sam --inputModel modelD2-HQ2D-1.hmm --jobTree ./jobTree-h1
	marginAlign/marginAlign D2/HQ2D/HQ2DReads_2.fasta D5/U00096.3.fasta outD2-HQ2D-100.sam --inputModel modelD2-HQ2D-100.hmm --jobTree ./jobTree-h100

	marginAlign/marginAlign D2/Normal2D/Nr2DReads_2.fasta D5/U00096.3.fasta outD2-Nr2D-1.sam --inputModel modelD2-Nr2D-1.hmm --jobTree ./jobTree-n1
	marginAlign/marginAlign D2/Normal2D/Nr2DReads_2.fasta D5/U00096.3.fasta outD2-Nr2D-100.sam --inputModel modelD2-Nr2D-100.hmm --jobTree ./jobTree-n100

	marginAlign/marginAlign D2/Reverse/ReverseReads_2.fasta D5/U00096.3.fasta outD2-Rev-1.sam --inputModel modelD2-Rev-1.hmm --jobTree ./jobTree-r1
	marginAlign/marginAlign D2/Reverse/ReverseReads_2.fasta D5/U00096.3.fasta outD2-Rev-100.sam --inputModel modelD2-Rev-100.hmm --jobTree ./jobTree-r100


	# D3 train: 
	#-----------
	# D3-609, 1%
	marginAlign/marginAlign D3/ERX955609/ERX955609_1_1_group.fastq D5/U00096.3.fasta paramsD3-609-1.sam --em --outputModel modelD3-609-1.hmm --jobTree ./jobTree609-1
	# D3-609, 100%
	marginAlign/marginAlign D3/ERX955609/ERX955609_1.fastq D5/U00096.3.fasta paramsD3-609-100.sam --em --outputModel modelD3-609-100.hmm --jobTree ./jobTree609-100
	
	# D3-580, 1%
	marginAlign/marginAlign D3/ERX963580/ERX963580_1_1_group.fastq D5/U00096.3.fasta paramsD3-580-1.sam --em --outputModel modelD3-580-1.hmm --jobTree ./jobTree580-1
	# D3-580, 100%
	marginAlign/marginAlign D3/ERX963580/ERX963580_1.fastq D5/U00096.3.fasta paramsD3-580-100.sam --em --outputModel modelD3-580-100.hmm --jobTree ./jobTree580-100

	# D3-855, 1%
	marginAlign/marginAlign D3/ERX978855/ERX978855_1_1_group.fastq D5/U00096.3.fasta paramsD3-855-1.sam --em --outputModel modelD3-855-1.hmm --jobTree ./jobTree855-1
	# D3-855, 100%
	marginAlign/marginAlign D3/ERX978855/ERX978855_1.fastq D5/U00096.3.fasta paramsD3-855-100.sam --em --outputModel modelD3-855-100.hmm --jobTree ./jobTree855-100


	# D3 align: 
	#-----------     
	marginAlign/marginAlign D3/ERX955609/ERX955609_2.fastq D5/U00096.3.fasta outD3-609-1.sam --inputModel modelD3-609-1.hmm --jobTree ./jobTree609-1
	marginAlign/marginAlign D3/ERX955609/ERX955609_2.fastq D5/U00096.3.fasta outD3-609-100.sam --inputModel modelD3-609-100.hmm --jobTree ./jobTree609-100

	marginAlign/marginAlign D3/ERX963580/ERX963580_2.fastq D5/U00096.3.fasta outD3-580-1.sam --inputModel modelD3-580-1.hmm --jobTree ./jobTree580-1
	marginAlign/marginAlign D3/ERX963580/ERX963580_2.fastq D5/U00096.3.fasta outD3-580-100.sam --inputModel modelD3-580-100.hmm --jobTree ./jobTree580-100

	marginAlign/marginAlign D3/ERX978855/ERX978855_2.fastq D5/U00096.3.fasta outD3-855-1.sam --inputModel modelD3-855-1.hmm --jobTree ./jobTree855-1
	marginAlign/marginAlign D3/ERX978855/ERX978855_2.fastq D5/U00096.3.fasta outD3-855-100.sam --inputModel modelD3-855-100.hmm --jobTree ./jobTree855-100


	# D4 train: 
	#-----------
	# D4-MAP, 1%
	marginAlign/marginAlign D4/MAP006-1_1_1_group.fasta D5/U00096.3.fasta paramsD4-MAP-1.sam --em --outputModel modelD4-MAP-1.hmm --jobTree ./jobTree-m1
	# D4-MAP, 100% 
	marginAlign/marginAlign D4/MAP006-1_1.fasta D5/U00096.3.fasta paramsD4-MAP-100.sam --em --outputModel modelD4-MAP-100.hmm --jobTree ./jobTree-m100


	# D4 align: 
	#-----------     
	marginAlign/marginAlign D4/MAP006-1_2.fasta D5/U00096.3.fasta outD4-MAP-1.sam --inputModel modelD4-MAP-1.hmm --jobTree ./jobTree-m1
	marginAlign/marginAlign D4/MAP006-1_2.fasta D5/U00096.3.fasta outD4-MAP-100.sam --inputModel modelD4-MAP-100.hmm --jobTree ./jobTree-m100


# Step 8:
# Repeat step 6 for LAST 
last_train_align:
	# Syntax:
	# lastdb -cR01 humdb humanMito.fa
	# last-train mydb queries.fasta
	# lastal humdb fuguMito.fa > myalns.maf

	# create db
	lastdb -cR01 ecolidb D5/U00096.3.fasta 

	#d1_train&align:
	# 1%
	last-train ecolidb D1/ERX708228/ERX708228_1_1_group.fastq	
	lastal ecolidb D1/ERX708228/ERX708228_2.fastq > alignD1-228-1.maf

	last-train ecolidb D1/ERX708229/ERX708229_1_1_group.fastq
	lastal ecolidb D1/ERX708229/ERX708229_2.fastq > alignD1-229-1.maf
	
	last-train ecolidb D1/ERX708229/ERX708230_1_1_group.fastq
	lastal ecolidb D1/ERX708230/ERX708230_2.fastq > alignD1-230-1.maf
	
	last-train ecolidb D1/ERX708229/ERX708231_1_1_group.fastq
	lastal ecolidb D1/ERX708231/ERX708231_2.fastq > alignD1-231-1.maf

	# 100%
	last-train ecolidb D1/ERX708228/ERX708228_1.fastq	
	lastal ecolidb D1/ERX708228/ERX708228_2.fastq > alignD1-228-100.maf

	last-train ecolidb D1/ERX708229/ERX708229_1.fastq
	lastal ecolidb D1/ERX708229/ERX708229_2.fastq > alignD1-229-100.maf
	
	last-train ecolidb D1/ERX708229/ERX708230_1.fastq
	lastal ecolidb D1/ERX708230/ERX708230_2.fastq > alignD1-230-100.maf
	
	last-train ecolidb D1/ERX708229/ERX708231_1.fastq
	lastal ecolidb D1/ERX708231/ERX708231_2.fastq > alignD1-231-100.maf

	#d2_train&align:
	# 1%
	last-train ecolidb D2/Forward/ForwardReads_1_1_group.fasta
	lastal ecolidb D2/Forward/ForwardReads_2.fasta > alignD2-For-1.maf

	last-train ecolidb D2/HQ2D/HQ2DReads_1_1_group.fasta
	lastal ecolidb D2/HQ2D/HQ2DReads_2.fasta > alignD2-HQ2D-1.maf
	
	last-train ecolidb D2/Normal2D/Nr2DReads_1_1_group.fasta
	lastal ecolidb D2/Normal2D/Nr2DReads_2.fasta > alignD2-Nr2D-1.maf
	
	last-train ecolidb D2/Reverse/ReverseReads_1_1_group.fasta
	lastal ecolidb D2/Reverse/ReverseReads_2.fasta > alignD2-Rev-1.maf

	# 100%
	last-train ecolidb D2/Forward/ForwardReads_1.fasta
	lastal ecolidb D2/Forward/ForwardReads_2.fasta > alignD2-For-100.maf

	last-train ecolidb D2/HQ2D/HQ2DReads_1.fasta
	lastal ecolidb D2/HQ2D/HQ2DReads_2.fasta > alignD2-HQ2D-100.maf
	
	last-train ecolidb D2/Normal2D/Nr2DReads_1.fasta
	lastal ecolidb D2/Normal2D/Nr2DReads_2.fasta > alignD2-Nr2D-100.maf
	
	last-train ecolidb D2/Reverse/ReverseReads_1.fasta
	lastal ecolidb D2/Reverse/ReverseReads_2.fasta > alignD2-Rev-100.maf

	#d3_train&align:
	# 1%
	last-train ecolidb D3/ERX955609/ERX955609_1_1_group.fastq
	lastal ecolidb D3/ERX955609/ERX955609_2.fastq > alignD3-609-1.maf
	
	last-train ecolidb D3/ERX963580/ERX963580_1_1_group.fastq
	lastal ecolidb D3/ERX963580/ERX963580_2.fastq > alignD3-580-1.maf
	
	last-train ecolidb D3/ERX978855/ERX978855_1_1_group.fastq 
	lastal ecolidb D3/ERX978855/ERX978855_2.fastq > alignD3-855-1.maf

	# 100%
	last-train ecolidb D3/ERX955609/ERX955609_1.fastq
	lastal ecolidb D3/ERX955609/ERX955609_2.fastq > alignD3-609-100.maf
	
	last-train ecolidb D3/ERX963580/ERX963580_1.fastq
	lastal ecolidb D3/ERX963580/ERX963580_2.fastq > alignD3-580-100.maf
	
	last-train ecolidb D3/ERX978855/ERX978855_1.fastq 
	lastal ecolidb D3/ERX978855/ERX978855_2.fastq > alignD3-855-100.maf

	#d4_train&align:
	# 1%
	last-train ecolidb D4/MAP006-1_1_1_group.fasta
	lastal ecolidb D4/MAP006-1_2.fasta > alignD4-MAP-1.maf

	# 100%
	last-train ecolidb D4/MAP006-1_1.fasta
	lastal ecolidb D4/MAP006-1_2.fasta > alignD4-MAP-100.maf


# Step 9:
# Repeat step 6 for BWA 
bwa_align:
	# Syntax:
	# bwa index ref.fa
	# bwa mem ref.fa reads.fq > aln-se.sam

	# create index
	bwa-0.7.12/bwa index D5/U00096.3.fasta

	#d1_align:
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D1/ERX708228/ERX708228_2.fastq > alignD1-228.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D1/ERX708229/ERX708229_2.fastq > alignD1-229.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D1/ERX708230/ERX708230_2.fastq > alignD1-230.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D1/ERX708231/ERX708231_2.fastq > alignD1-231.sam

	#d2_align:	
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D2/Forward/ForwardReads_2.fasta > alignD2-For.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D2/HQ2D/HQ2DReads_2.fasta > alignD2-HQ2D.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D2/Normal2D/Nr2DReads_2.fasta > alignD2-Nr2D.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D2/Reverse/ReverseReads_2.fasta > alignD2-Rev.sam

	#d3_align:	
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D3/ERX955609/ERX955609_2.fastq > alignD3-609.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D3/ERX963580/ERX963580_2.fastq > alignD3-580.sam
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D3/ERX978855/ERX978855_2.fastq > alignD3-855.sam

	#d4_align:	
	bwa-0.7.12/bwa mem D5/U00096.3.fasta D4/MAP006-1_2.fasta > alignD4-MAP.sam


