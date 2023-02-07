# -*- coding: utf-8 -*-

#Created on January 2023

#@author: Xiaokun Hong,Jiyang Lv


'''
main program for hosting the website PhosPPI
'''

from flask import Flask, render_template, request, send_from_directory
from flask_mail import Mail,Message
from flask_bootstrap import Bootstrap
from flask_wtf import FlaskForm
from wtforms import SelectField, SubmitField, RadioField, StringField, BooleanField,FileField
from wtforms.validators import DataRequired
import os
import pandas as pd
import zipfile
from flask import jsonify
from flask import jsonify,make_response
import json
from werkzeug.utils import secure_filename
import numpy as np
import time
import shutil

import pandas as pd
import numpy as np
import sys
import sklearn
from sklearn.metrics import accuracy_score,precision_score,recall_score,f1_score,roc_auc_score
from sklearn.metrics import confusion_matrix
import time
import pickle
import lightgbm as lgb
import traceback

import datetime

import random
import string

import sqlite3
from threading import Thread
from flask_executor import Executor





MODEL1= pickle.load(open("model1.dat","rb"))
MODEL2= pickle.load(open("model2.dat","rb"))

COUNT = int(open('static/count.txt').read())
job_number=24

max_running_number=5
max_user_queued_number=9

#####################网页相关####################################################



conn = sqlite3.connect('job3.db',check_same_thread=False)
c = conn.cursor()
#print ("数据库打开成功")
c.execute("""CREATE TABLE IF NOT EXISTS mytable3 (iden TEXT PRIMARY KEY,
    status TEXT,
    submitted_time TEXT,
    job_name TEXT, 
    job_number2 TEXT, 
    sitenumber INT, 
    running_time TEXT, 
    ww_list TEXT, 
    A1 TEXT, 
    B1 TEXT, 
    C1 TEXT, 
    D1 TEXT, 
    E1 TEXT, 
    F1 TEXT, 
    zipname TEXT, 
    mm2 TEXT, 
    nn INT, 
    ns INT, 
    nt INT, 
    ny INT,
    email TEXT,
    SITE TEXT,
    start REAL)""")
#c.execute("alter table mytable3 add SITE text")
conn.commit()


class Config:
    SECRET_KEY = 'hard to guess string'
    SSL_DISABLE = False
    WTF_CSRF_ENABLED = False
    DEBUG = False


class NonValidatingSelectField(SelectField):
    # Skip the pre-validation, otherwise it will raise the "Not a valid choice" error
    def pre_validate(self, form):
        pass



class QForm(FlaskForm):
    emails = StringField()
    jobs = StringField()
    P1 = StringField()
    P2 = StringField()
    SITE = StringField()
   

    
    
    submit = SubmitField('Search')


#####################主函数####################################################
app = Flask(__name__)
app.config.from_object(Config)
app.config['MAIL_SERVER']='smtp.qq.com'         #邮件服务器的名称/IP地址
app.config['MAIL_PORT'] = 465                  #所用服务器的端口号
app.config['MAIL_USERNAME'] = 'hxk33461@foxmail.com'     #发件人的用户名
app.config['MAIL_PASSWORD'] = 'fhdoocpoyvrubchf'         #发件人的POP3/IMAP/SMTP服务的SSL连接客户端授权码
app.config['MAIL_USE_TLS'] = False              #禁用传输安全层加密
app.config['MAIL_USE_SSL'] = True               #启用安全套接字层加密



app.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024 #只能上传最大2M的文件
#app.config['UPLOAD_FOLDER']= './upload'

executor = Executor(app)

mail = Mail(app)                                #创建邮件类对象


Bootstrap(app)

def global_count():
    global COUNT
    return COUNT



#产生唯一标识符
def random_string_generator(str_size, allowed_chars):
    return ''.join(random.choice(allowed_chars) for x in range(str_size))
    
    

    
app.add_template_global(global_count,'global_count')

#判断序列是否为标准氨基酸
def sequences_extract(chainA):
    target_sequence1=''
    target_sequence2=''
    target_sequence3=''
    target_sequence4=''
        
    for line in chainA:
        target_sequence1 = target_sequence1 + line
    
    #去掉FASTA文件标题行，并且将每行的序列拼接起来
    target_sequence1=target_sequence1.split('\n')
    
    for k in range(1,len(target_sequence1)):
        target_sequence2=target_sequence2+target_sequence1[k].strip()
       
  
    #检查是否是标准氨基酸，并将全部氨基酸转换为大写
    for j in target_sequence2:
        if j.upper() not in 'ARNDCQEGHILKMFPSTWYV':
            target_sequence3=''
            break
        else:
           target_sequence3=target_sequence3+j.upper()
    target_sequence4=target_sequence1[0]+"\n"+target_sequence3
    return target_sequence3,target_sequence4
    



 
    
#实现对用户上传fasta文件的处理
def fasta_sequences_extract1(PATH,FILE, resid):
    resid = int(resid) - 1
    target_sequence = ''
    sequence=''
    
    paths=os.path.join(PATH, FILE)
    for line in open(paths):
        str = line.strip()
        if str[0] != '>':
            target_sequence = target_sequence + str
    for i in range(int(resid) - 15, int(resid) + 16):
        if i < 0:
            sequence=sequence+"*"
        else:
           try:
              target_sequence[i]
           except IndexError:
              sequence=sequence+"*"
           else:
              sequence=sequence+target_sequence[i]
    return sequence


def fasta_sequences_extract2(PATH,FILE, resid):
    resid = int(resid) - 1
    target_sequence = ''
    sequence=''
    
    paths=os.path.join(PATH, FILE)
    for line in open(paths):
        str = line.strip()
        if str[0] != '>':
            target_sequence = target_sequence + str
    for i in range(int(resid) - 5, int(resid) + 4):
        if i < 0:
            sequence=sequence+"*"
        else:
           try:
              target_sequence[i]
           except IndexError:
              sequence=sequence+"*"
           else:
              sequence=sequence+target_sequence[i]
    return sequence



def check_site(PATH,FILE,resid):
    resid = int(resid) - 1
    target_sequence = ''
    sequence=''
    
    paths=os.path.join(PATH, FILE)
    for line in open(paths):
        str = line.strip()
        if str[0] != '>':
            target_sequence = target_sequence + str
            
    amino=target_sequence[resid]
             
    return amino 
    

def sequence_all(PATH,FILE):
    target_sequence = ''  
    paths=os.path.join(PATH, FILE)
    for line in open(paths):
        str = line.strip()
        if str[0] != '>':
            target_sequence = target_sequence + str         
    return target_sequence 


  
#定义氨基酸标签OBC 
def aa_label(three_letter):
    aa = {'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  'R':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  'N':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  'D':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  'C':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  
    'Q':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  'E':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  'G':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],  'H':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],  'I':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],  
    'L':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],  'K':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],  'M':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],  'F':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],  'P':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
    'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],  'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],  'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],  'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],  'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],  '*':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]}
    return aa.get(three_letter, 0)    
    

#para1-#CIDH920105 Normalized average hydrophobicity scales (Cid et al., 1992)
def parameter1(three_letter):
    aa = {'A':0.02,  'R':-0.42,  'N':-0.77,  'D':-1.04,  'C':0.77,  
    'Q':-1.10,  'E':-1.14,  'G':-0.80,  'H':0.26,  'I':1.81,  
    'L':1.14,  'K':-0.41,  'M':1.00,  'F':1.35,  'P':-0.09,
    'S':-0.97,  'T':-0.77,  'W':1.71,  'Y':1.11,  'V':1.13,  '*':0}
    return aa.get(three_letter, 0)

#para2-#BHAR880101 Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)
def parameter2(three_letter):
    aa = {'A':0.357,  'R':0.529,  'N':0.463,  'D':0.511,  'C':0.346,  
    'Q':0.493,  'E':0.497,  'G':0.544,  'H':0.323,  'I':0.462,  
    'L':0.365,  'K':0.466,  'M':0.295,  'F':0.314,  'P':0.509,
    'S':0.507,  'T':0.444,  'W':0.305,  'Y':0.420,  'V':0.386,  '*':0}
    return aa.get(three_letter, 0)

#para3-#CHAM820101 Polarizability parameter (Charton-Charton, 1982)
def parameter3(three_letter):
    aa = {'A':0.046,  'R':0.291,  'N':0.134,  'D':0.105,  'C':0.128,  
    'Q':0.180,  'E':0.151,  'G':0.000,  'H':0.230,  'I':0.186,  
    'L':0.186,  'K':0.219,  'M':0.221,  'F':0.290,  'P':0.131,
    'S':0.062,  'T':0.108,  'W':0.409,  'Y':0.298,  'V':0.140,  '*':0}
    return aa.get(three_letter, 0)

#para4-#CHAM820102 Free energy of solution in water, kcal/mole (Charton-Charton, 1982)
def parameter4(three_letter):
    aa = {'A':-0.368,  'R':-1.03,  'N':0.,  'D':2.06,  'C':4.53,  
    'Q':0.731,  'E':1.77,  'G':-0.525,  'H':0.,  'I':0.791,  
    'L':1.07,  'K':0.,  'M':0.656,  'F':1.06,  'P':-2.24,
    'S':-0.524,  'T':0.,  'W':1.60,  'Y':4.91,  'V':0.401,  '*':0}
    return aa.get(three_letter, 0)

#para5-#CHOC760101 Residue accessible surface area in tripeptide (Chothia, 1976)
def parameter5(three_letter):
    aa = {'A':115.,  'R':225.,  'N':160.,  'D':150.,  'C':135.,  
    'Q':180.,  'E':190.,  'G':75.,  'H':195.,  'I':175.,  
    'L':170.,  'K':200.,  'M':185.,  'F':210.,  'P':145.,
    'S':115.,  'T':140.,  'W':255.,  'Y':230.,  'V':155.,  '*':0}
    return aa.get(three_letter, 0)

#para6-BIGC670101 Residue volume (Bigelow, 1967)
def parameter6(three_letter):
    aa = {'A':52.6,  'R':109.1,  'N':75.7,  'D':68.4,  'C':68.3,  
    'Q':89.7,  'E':84.7,  'G':36.3,  'H':91.9,  'I':102.0,  
    'L':102.0,  'K':105.1,  'M':97.7,  'F':113.9,  'P':73.6,
    'S':54.9,  'T':71.2,  'W':135.4,  'Y':116.2,  'V':85.1,  '*':0}
    return aa.get(three_letter, 0)
      
#para7-#CHAM810101 Steric parameter (Charton, 1981)
def parameter7(three_letter):
    aa = {'A':0.52,  'R':0.68,  'N':0.76,  'D':0.76,  'C':0.62,  
    'Q':0.68,  'E':0.68,  'G':0.00,  'H':0.70,  'I':1.02,  
    'L':0.98,  'K':0.68,  'M':0.78,  'F':0.70,  'P':0.36,
    'S':0.53,  'T':0.50,  'W':0.70,  'Y':0.70,  'V':0.76,  '*':0}
    return aa.get(three_letter, 0)
    
#para8-#DAYM780201 Relative mutability (Dayhoff et al., 1978b)
def parameter8(three_letter):
    aa = {'A':100.,  'R':65.,  'N':134.,  'D':106.,  'C':20.,  
    'Q':93.,  'E':102.,  'G':49.,  'H':66.,  'I':96.,  
    'L':40.,  'K':56.,  'M':94.,  'F':41.,  'P':56.,
    'S':120.,  'T':97.,  'W':18.,  'Y':41.,  'V':74.,  '*':0}
    return aa.get(three_letter, 0)
  
#para9-#GRAR740102 Polarity (Grantham, 1974)
def parameter9(three_letter):
    aa = {'A':8.1,  'R':10.5,  'N':11.6,  'D':13.0,  'C':5.5,  
    'Q':10.5,  'E':12.3,  'G':9.0,  'H':10.4,  'I':5.2,  
    'L':4.9,  'K':11.3,  'M':5.7,  'F':5.2,  'P':8.0,
    'S':9.2,  'T':8.6,  'W':5.4,  'Y':6.2,  'V':5.9,  '*':0}
    return aa.get(three_letter, 0)

#2>/dev/null
def submit_pssm(PATH2,file_name):

    os.system('/array/zxli/HXK_filefolder/PhosPPI/ncbi-blast-2.9.0+/bin/psiblast -db '+'/array/zxli/HXK_filefolder/PhosPPI/zypnr/nr '+'-num_iterations 3 -evalue 0.001 -query '+PATH2+file_name+' -out_ascii_pssm '+PATH2+file_name + '.pssm >>TMP.txt')



#实现对psi-blast产生的PSSM文件的处理
def readPssm(PATH,pssmFile):
	# index of 'ACDE..' in 'ARNDCQEGHILKMFPSTWYV'(blast order)
	paths=os.path.join(PATH, pssmFile)
	fp = open(paths, 'r')
	lines = fp.readlines()
	fp.close()
	pssm = []
	#print(lines)	
	for line in lines:
		splitLine = line.split()
		# valid lines should have 32 points of data.
		# any line starting with a # is ignored
		if (len(splitLine) == 44) and (splitLine[0] != '#'):
			pssmTemp = [float(1.0/(1.0 + np.exp(-float(i)))) for i in splitLine[2:22]]
			pssm.append(pssmTemp)
	pssm1 = np.array(pssm)
	mean_pssm = pssm1.mean(axis=0)
	mean_pssm1 = list(mean_pssm)
	return mean_pssm1

#计算netsurfp特征
def submit_netsurfp(PATH2,PATH,file_name):   
    os.system("python ./NetSurfP-3.0_standalone/nsp3.py -m ./NetSurfP-3.0_standalone/models/nsp3.pth -i "+PATH2+file_name+" -o "+PATH+" >>TMP.txt")
    
#拷贝netsurfp结果文件    
def get_dir(PATH,path,fileType,job_number2):
    allfilelist=os.listdir(path)
    
    for file in allfilelist:
        filepath =os.path.join(path,file)
        
        if os.path.isdir(filepath):
            allfilelist2 = os.listdir(filepath)
        
            for file2 in allfilelist2:
                filepath3 = os.path.join(filepath, file2)
               
                if filepath3.endswith(fileType):
                    #print(filepath3)
                    shutil.copy(filepath3, PATH)                                
        else:
            pass
            #print("No csv file")
           
           
    
    job_path1=PATH+"/*.csv"
    job_path2=PATH+"/"+job_number2+".csv"
    job_path3=PATH+"/01"

    os.system("mv %s %s"%(job_path1,job_path2))

    os.system("rm -f -r %s"%job_path3)
    os.system("rm -f TMP.txt")



#提取netsurfp特征-ASA,SS
def extract_netsurfp(filename,resid):
    resid = int(resid) - 1
    ASA=[]
    SS=[]
    df=pd.read_csv(filename,sep=",")
    #提取ASA
    for i in range(resid-15,resid+16):
       if i <0:
          ASA.append(0)
       else:
          try:
             df.iloc[i,4]
          except IndexError:
             ASA.append(0)
          else:
             ASA.append(df.iloc[i,4])
    #提取SS      
    for i in range(resid-15,resid+16):
       if i <0:
          SS.append(0)
          SS.append(0)
          SS.append(0)   
       else:
         try:
             df.iloc[i,6]
         except IndexError:
             SS.append(0)
             SS.append(0)
             SS.append(0)
         else:
             SS.append(df.iloc[i,6])
             SS.append(df.iloc[i,7])
             SS.append(df.iloc[i,8]) 
       
    return ASA,SS    
     


#定义上传的文件类型
ALLOWED_EXTENSIONS = set(['fasta','FASTA'])
#UPLOAD_FOLDER = './upload'

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS
           
           


# Extra redirect request for "/favicon.ico"
@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.before_request
def count():
    global COUNT
    import re
    if not re.match(r'^.*\.\w{1,6}$',request.path):
        COUNT += 1
    if COUNT % 10 ==0:
        open('static/count.txt','w').write(str(COUNT))
        

    
# Home page
@app.route('/')
def index():
    return render_template('index.html')

# Other pages
@app.route('/<string:page>')
def show(page):
    return render_template('{}.html'.format(page))
    
# data
@app.route('/data/<string:page>')
def data(page):
    return send_from_directory('data/', page)




@app.route('/Home', methods=['GET', 'POST'])

def index2():

    global job_number,MODEL1,MODEL2,db,Q,mail
    if request.method == 'POST':
        job_number=job_number+1
        job_number2 =str('%08d' %job_number)
        
        PATH="./"+job_number2
        PATH2="./"+job_number2+"/"
        
        #获取表单数据
        form = QForm(data=request.form)
        
        P1=form.P1.data
        P2=form.P2.data
        
        SITE=form.SITE.data
        site_list = SITE.split(',')
        
        emailname = form.emails.data
        
        job_name = form.jobs.data
           
        #判断是否为标准氨基酸
        SEQ1,SEQ11=sequences_extract(P1)
        SEQ2,SEQ22=sequences_extract(P2)
        
        #        
        file1 = request.files['file1']
        #filename1 = secure_filename(file1.filename)
        filename1="P1.fasta" 
        
        #file2 = request.files['file2']
        #filename2 = secure_filename(file2.filename)
                   
        file3 = request.files['file3']
        #filename3 = secure_filename(file3.filename)
        filename3="P2.fasta" 
        
        
        
        #if not (email and theme and SITE and bool(SEQ1) ^ bool(file1.filename) and bool(SEQ2) ^ bool(file3.filename)):
            #return render_template('index.html',XX="Please complete all fields correctly!")

        if not (emailname):
            return render_template('index.html',XX="Please input email address!")
       
        if not (job_name):
            return render_template('index.html',XX="Please input job name!")

        if not (SITE):
            return render_template('index.html',XX="Please designate positions of phosphorylation sites!")
  
        if not (bool(SEQ1) ^ bool(file1.filename)):
            return render_template('index.html',XX="Please only paste P1 protein sequence or upload FASTA file!")
    
        if not (bool(SEQ2) ^ bool(file3.filename)):
            return render_template('index.html',XX="Please only paste P2 protein sequence or upload FASTA file!")
             
                                          
        os.system("mkdir %s " %job_number2)

       
       #保存文件
        if file1.filename:  
            if allowed_file(file1.filename):
                file1.save(os.path.join(PATH, filename1))    
            else:
                return render_template('index.html',XX="Please input FASTA file correctly!") 
        else:
            paths1=os.path.join(PATH,filename1)      
            open(paths1,'w',newline="").write((P1+"\n") if P1[-1]!="\n" else P1) #判断最后一个字符是否为分行符号，如无则增加一个分行符号

        if file3.filename:
            if allowed_file(file3.filename):
                file3.save(os.path.join(PATH, filename3))
            else:
                return render_template('index.html',XX="Please input FASTA file correctly!") 
        else:
            paths3=os.path.join(PATH,filename3)      
            open(paths3,'w',newline="").write((P2+"\n") if P2[-1]!="\n" else P2) #判断最后一个字符是否为分行符号，如无则增加一个分行符号                          
        
        for site in SITE.split(','):
            amino=check_site(PATH,filename1,site)
            if amino not in ['S','T','Y']:
                return render_template('index.html',XX="Please designate positions of phosphorylation sites correctly!")
        
        curr_time = datetime.datetime.now()
        start=time.time()
        submitted_time = datetime.datetime.strftime(curr_time,'%Y-%m-%d %H:%M:%S')
        identifier=random_string_generator(10, string.ascii_letters + string.digits)+str(int(start))
        #数据库增加数据

        c.execute("select count(1) from mytable3 where email = ? and status in ('waiting','running')",(emailname,))
        r = c.fetchone()
        if r and r[0]>max_user_queued_number:
            return render_template('index.html',XX="your queued task number is greater than 10!") 
            
        c.execute("INSERT INTO mytable3 VALUES(?,'waiting',?,?,?,1,'1','1','1','1','1','1','1','1','1','1',1,1,1,1,?,?,?)",(identifier,submitted_time,job_name,job_number2,emailname,SITE,start))   
        conn.commit()
        
        theme="PhosPPI Server results for job "+ job_name
        Diag="Dear researcher:<br/>Congratulations!<br/>Your job has been submitted successfully!<br/>Thank you for using PhosPPI Server!<br/>PhosPPI Server"
        msg = Message(theme, sender='hxk33461@foxmail.com', recipients=[emailname],html=Diag)
        mail.send(msg) 
        
        c.execute("select count(1) from mytable3 where status in ('running')")
        r = c.fetchone()
        if not r or r[0]<max_running_number:
            executor.submit(create_task,identifier)
                                       
        #return render_template('result3.html',sitenumber=sitenumber,running_time=running_time,ww_list=ww_list,A1=A1,B1=B1,C1=C1,D1=D1,E1=E1,F1=F1,zipname=zipname,mm2=mm2,nn=nn,ns=ns,nt=nt,ny=ny)
    return render_template('index.html')

def create_task(identifier):
    conn = sqlite3.connect('job3.db',check_same_thread=False)
    c = conn.cursor()
    c.execute("update mytable3 set status='running' where iden = ?",(identifier,))   
    conn.commit()
    try:
        c.execute("select * from mytable3 where iden = ?",(identifier,))
        r = c.fetchone()
        
        PSSM1, PSSM3 = [], []
        
        job_number2=r[4]
        PATH="./"+job_number2
        PATH2="./"+job_number2+"/"    
        
        filename1="P1.fasta"
        filename3="P2.fasta"
        
        SITE=r[21]
        emailname=r[20]
        
        job_name=r[3]
        theme="PhosPPI Server results for job "+ job_name
        
        submitted_time=r[2]
        
        start=r[22]
        
        ww_list=[]
     
        
        for site in SITE.split(','):
            amino=check_site(PATH,filename1,site)                  
            array13, tt13, AAC, para1, para2, para3, para4, para5, para6, para7, para8, para9 = [],[],[],[],[],[],[],[],[],[],[],[]
      
      #for model1 feature
            SEQ13=fasta_sequences_extract1(PATH,filename1,site)
                                           
            for i in SEQ13:
                array13.append(i)
            for i in range(len(array13)):
                tt13.append(aa_label(array13[i]))
          
                para1.append(parameter1(array13[i]))
                para2.append(parameter2(array13[i]))
                para3.append(parameter3(array13[i]))
                para4.append(parameter4(array13[i]))
                para5.append(parameter5(array13[i]))
                para6.append(parameter6(array13[i]))
                para7.append(parameter7(array13[i]))
                para8.append(parameter8(array13[i]))
                para9.append(parameter9(array13[i]))
                
            AAC = [array13.count(i)/len(array13) for i in "ARNDCQEGHILKMFPSTWYV*"]
      
          #将OBC每个氨基酸的OBC合并    
            OBC13=[y for x in tt13 for y in x]

      #计算SS,ASA，并保存在新建文件夹中
            path11=os.path.join(PATH, "{}.csv".format(job_number2))
            if os.path.exists(path11):
                #path11=os.path.join(PATH, "{}.csv".format(job_number2))
                ASA,SS=extract_netsurfp(path11,site)                
            
            else:
                submit_netsurfp(PATH2,PATH,filename1)  
                PATH99=PATH2+"01"
                get_dir(PATH,PATH99,".csv",job_number2)                              
                #path11=os.path.join(PATH, "{}.csv".format(job_number2))
                ASA,SS=extract_netsurfp(path11,site)
        
            Feature1=AAC+OBC13+para1+para2+para3+para4+para5+para6+para7+para8+para9+SS+ASA
            
            Feature1=pd.DataFrame(Feature1).T
            
            XX1=Feature1.iloc[:,sorted([1060,1059,1053,1054,15,1064,1055,1062,1073,1044,1048,1058,1068,371,1066,1065,1056,1061,1050,1052,1063,1074,1047,1046,1072,1070,1057,1071,354,748,905,781,1045,1000,778,274,1067,1043,10,746,1069,747,1051,685,18,1042,715,1049,309,903,952,1041,330,6,1,719,902,978,912,1036,723,762,906,936,1039,3,724,901,919,1015,711,868,932,951,964,1011,1027,331,981,985,954,1034,742,752,779,783,993,1024])].values
            Pred1=MODEL1.predict(XX1)#预测测试集分类
            Prob1=MODEL1.predict_proba(XX1)[:, 1] 
            
            
            if Pred1==0:
                Prob1=1-Prob1 
                Pred="No effect"
                Prob=np.round(Prob1,2)  
                
            else:
                array14=[]
                tt14=[]
              #model2-P1-pssm
                SEQ14=fasta_sequences_extract2(PATH,filename1,site)
                for i in SEQ14:
                    array14.append(i)                   
                for i in range(len(array14)):
                    tt14.append(aa_label(array14[i]))     
              #将OBC每个氨基酸的OBC合并  model2-P1-OBC  
                OBC14=[y for x in tt14 for y in x] 
                    
                if not PSSM1:
                    submit_pssm(PATH2,filename1)
                    PSSM1=readPssm(PATH,filename1+".pssm")        

                if not PSSM3:
                    
                    if sequence_all(PATH,filename1)==sequence_all(PATH,filename3):
                        PSSM3=PSSM1
                    else:
                        submit_pssm(PATH2,filename3)
                        PSSM3=readPssm(PATH,filename3+".pssm")

                Feature2=OBC14+PSSM1+PSSM3
                
                Feature2=pd.DataFrame(Feature2).T
                
                XX2=Feature2.iloc[:,:].values
                Pred2=MODEL2.predict(XX2)#预测测试集分类
                Prob2=MODEL2.predict_proba(XX2)[:, 1]
              
                if Pred2==0:
                    Prob2=1-Prob2 
                    Pred="Inhibit"
                    Prob=np.round(Prob2,2)   
                    
                else:
                    Pred="Enhance"   
                    Prob=np.round(Prob2,2)                    
                
            ww_list.append({'Position':site,'AA':amino,'Effect':Pred,'Score':Prob[0]})

       
       
        data=pd.read_csv(path11,sep=",")
        data.columns=['id','seq','n','rsa','asa','q3','p[q3_H]','p[q3_E]','p[q3_C]','q8','p[q8_G]','p[q8_H]','p[q8_I]','p[q8_B]','p[q8_E]','p[q8_S]','p[q8_T]','p[q8_C]','phi','psi','disorder']
        data['label']=data['seq']+data['n'].map(str)
       
        data['asa']=round(data['asa'],2)
        data['p[q3_H]']=round(data['p[q3_H]'],2)
        data['p[q3_E]']=round(data['p[q3_E]'],2)
        data['p[q3_C]']=round(data['p[q3_C]'],2)
        
        
        A1=data.iloc[:,2].to_list()
        B1=data.iloc[:,6].to_list()
        C1=data.iloc[:,7].to_list()
        D1=data.iloc[:,8].to_list()
       
        E1=A1
       
        F1=data.iloc[:,4].to_list()
       
       
       # markpoint
        
        mark=data[data['n'].isin([int(i) for i in SITE.split(",")])]
        
        mark1 = mark.iloc[:,[21,2,4]]
        mark1.columns = ['value', 'xAxis', 'yAxis']
        mm2 = [{'value': m['value'], 'xAxis': m['xAxis'], 'yAxis': m['yAxis'],
        'itemStyle': {'color': "rgba(115, 159, 250, .8)"}} for u, m in mark1.iterrows()]
       
        nn=data[(data['seq']=='S') | (data['seq']=='T') | (data['seq']=='Y')].shape[0]
        ns=data[(data['seq']=='S')].shape[0]
        nt=data[(data['seq']=='T')].shape[0]
        ny=data[(data['seq']=='Y')].shape[0]
        #print(ww_list)
        ww = pd.DataFrame(ww_list)
        tablename=PATH2+"Effect.csv"
        ww.to_csv(tablename,index=False)
       
        shutil.make_archive(PATH, 'zip', job_number2)
       
        zipname=job_number2 + ".zip"
                            
        os.system("mv %s ./data"%job_number2)
        os.system("mv %s ./data"%zipname)
       

        sitenumber=len(SITE.split(','))
                   
        end=time.time()
        running_time=end-start   
        running_time=round(running_time,2)
        
        
        #数据库更新数据
        c.execute("update mytable3 set status='completed',submitted_time=?,job_name=?, job_number2=?, sitenumber=?, running_time=?, ww_list=?, A1=?, B1=?, C1=?, D1=?, E1=?, F1=?, zipname=?, mm2=?, nn=?, ns=?, nt=?, ny=?,email=?,SITE=?,start=? where iden = ?",(submitted_time,job_name,job_number2,sitenumber,running_time,json.dumps(ww_list),json.dumps(A1),json.dumps(B1),json.dumps(C1),json.dumps(D1),json.dumps(E1),json.dumps(F1),zipname,json.dumps(mm2),nn,ns,nt,ny,emailname,SITE,start,identifier))
        conn.commit() 
      
       
        #邮件服务
        
        with app.app_context():
            url = "https://phosppi.sjtu.edu.cn/result/"+identifier
            Diag="Dear researcher:<br/>Congratulations<br/>Your job has been completed!<br/>The results can be retrieved from the following URL:<a href='"+url+"'>"+url+"</a><br/>Thank you for using PhosPPI Server!<br/>PhosPPI Server"
            msg = Message(theme, sender='hxk33461@foxmail.com', recipients=[emailname],html=Diag)
            mail.send(msg) 
    except:
        print(traceback.format_exc())

@app.route('/result/<string:identifier>')
def results(identifier):
    
    #conn = sqlite3.connect('job.db')
    #c = conn.cursor()
    
    database=c.execute("SELECT * FROM mytable3 WHERE iden == '%s' "%identifier)

    for xxx in database:
  
        sitenumber=xxx[5]
        running_time=xxx[6]
        ww_list=json.loads(xxx[7])
    
        A1=json.loads(xxx[8])
        B1=json.loads(xxx[9])
        C1=json.loads(xxx[10])
        D1=json.loads(xxx[11])
        E1=json.loads(xxx[12])
        F1=json.loads(xxx[13])
    
        zipname=xxx[14]
        mm2=json.loads(xxx[15])
    
        nn=xxx[16]
        ns=xxx[17]
        nt=xxx[18]
        ny=xxx[19]
    

    
    
    return render_template('result4.html',sitenumber=sitenumber,running_time=running_time,ww_list=ww_list,A1=A1,B1=B1,C1=C1,D1=D1,E1=E1,F1=F1,zipname=zipname,mm2=mm2,nn=nn,ns=ns,nt=nt,ny=ny)

  

@app.route('/Queue')
def results2():

    query_list=[]
    
    #conn = sqlite3.connect('job.db')
    #c = conn.cursor()
 
    db=c.execute("SELECT * FROM mytable3;")
    
    for row in db:
    
        job_number2=row[4]
        job_name=row[3]
        status=row[1]
        submitted_time=row[2]   
        query_list.append({'job_number2':job_number2,'job_name':job_name,'status':status,'submitted_time':submitted_time})
 
 
    #c.execute("SELECT * FROM mytable3;")
    #query_list=c.fetchall()
       
    return render_template('Queue.html',query=query_list)


def loop():
    conn = sqlite3.connect('job3.db',check_same_thread=False)
    c = conn.cursor()
    while 1:
        try:
            c.execute("select count(1) from mytable3 where status in ('running')")
            r = c.fetchone()
            #print(r)
            if not r or r[0]<max_running_number:
                c.execute("select * from mytable3 where status in ('waiting') order by submitted_time limit 1")
                res = c.fetchone()
                if res:
                    create_task(res[0])
        except:
            print(traceback.format_exc())
        time.sleep(2)


if __name__ == '__main__':
    # 对应intenal IP
    for i in range(max_running_number):
        Thread(target=loop,args=()).start()    
    app.run(host='0.0.0.0', port=5000, debug=False)



