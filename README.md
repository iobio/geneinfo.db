## Introduction
### This project is for rebuiding gene.iobio database based on [gene.model.db](https://github.com/iobio/gene.model.db/tree/master)

### Since the gff files are massive, it will take a while to process and read into the database. That said, there is space to optimize. For example, use Dask to process data parallelly.

____

## Instruction
### 1. Install 
``` git clone https://github.com/iobio/geneinfo.db.git``` <br>
``` cd your project directory```
### 2. Create your virtualenv
``` python3 -m venv myenv ```
``` source myenv/bin/activate ```
### 3. Install packages and dependencies
``` pip3 install -r requirements.txt```
### 4. Run
``` Python3 data_downloader.py``` <br>
``` Python3 database_builder.py```
### 5. To create a new genes.json which is used in the gene.iobio search type-ahead, run the following:
``` python3 create_genes.json``` <br>
``` cp genes.json ../gene.iobio/client/data/```

