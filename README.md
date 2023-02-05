# PhosPPI
A web service (PhosPPI) for **Sequence-based Machine Learning Method for Predicting the Effects of Phosphorylation on Protein-Protein Interactions** created by *Xiaokun Hong, Jiyang Lv*

## Set up environment

1. Our Conda environment was created based on NetSurfP-3.0's Conda environment (Nucleic Acids Res,2022).
	*Download links(Linux 3.0):https://services.healthtech.dtu.dk/service.php?NetSurfP-3.0
		* conda env create --file environment.yml
		* conda activate nsp3
		* python setup.py install
		* pip3 install -r requirements.txt

In addition, Web service and Machine Learning Python packages are required to install according to the following instructions:

2. Install Web service Python packages
	* pip install 
		* Flask=2.0.2
		* Flask-WTF=1.0.0
		* WTForms=3.0.1
		* Flask-Bootstrap=3.3.7.1
		* Jinja2=3.0.3
		* pandas=1.4.1
		* gunicorn=20.1.0
		* plotnine=0.8.0
		
		
3. Other packages no need to install (*latest version up to 2022-02-16*; best to install them as network connection from China may fail)
   * [Bootstrap](https://getbootstrap.com/): v5.1.3
   * [jQuery](https://jquery.com/): v3.6.0
   * [jQuery UI](https://jqueryui.com/): v1.13.1

4. To perform PSI-BLAST to obtain PSSM features, ncbi-blast-2.9.0+-x64-linux are also needed to install.

5. Machine Learning Python packages
	* pip install	
		-

