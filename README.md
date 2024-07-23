# PhosPPI
A web service (PhosPPI) for **Sequence-based Machine Learning Method for Predicting the Effects of Phosphorylation on Protein-Protein Interactions** created by *Xiaokun Hong, Jiyang Lv*

Hong X, Lv J, Li Z, Xiong Y, Zhang J*, Chen HF*. Sequence-based machine learning method for predicting the effects of phosphorylation on protein-protein interactions. Int J Biol Macromol. 2023 Jul 15;243:125233. doi: 10.1016/j.ijbiomac.2023.125233. Epub 2023 Jun 6. PMID: 37290543.

## Set up environment

1. Our Conda environment was created based on NetSurfP-3.0's Conda environment (Nucleic Acids Res,2022).
*https://services.healthtech.dtu.dk/service.php?NetSurfP-3.0
  * conda env create --file environment.yml
  * conda activate nsp3
  * python setup.py install
  * pip3 install -r requirements.txt

*In addition, Web service and Machine Learning Python packages are required to install according to the following instructions:

2. Install Web Service Python packages
* pip install 
  * Flask=1.1.2
  * Flask-WTF=1.1.1
  * WTForms=3.0.1
  * Flask-Bootstrap=3.3.7.1
  * Jinja2=2.11.3
  * pandas=1.2.3
  * gunicorn=19.9.0
		
		
3. Other packages no need to install (*latest version up to 2022-02-16*; best to install them as network connection from China may fail)
  * [Bootstrap](https://getbootstrap.com/): v5.1.3
  * [jQuery](https://jquery.com/): v3.6.0
  * [jQuery UI](https://jqueryui.com/): v1.13.1

4. To perform PSI-BLAST to obtain PSSM features, ncbi-blast-2.9.0+-x64-linux are also needed to install.
*https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/

5. Install Machine Learning Python packages
* pip install sklearn=0.0
* pip install lightgbm=3.3.3



